/**
 * @file pulseqlib_bridge.c
 * @brief POSIX bridge to pypulseq_host child process.
 *
 * Spawns the child via fork/exec with stdin/stdout pipes.
 * Communicates using the NimPulseqGUI persistent text protocol:
 *   LIST_PROTOCOL, VALIDATE, GENERATE, QUIT.
 *
 * Pure C89 + POSIX (no vendor dependencies).
 */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200112L
#endif

#include "pulseqlib_bridge.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

/* Avoid <sys/wait.h> and <signal.h> — the GE host32 cross-compiler
 * sysroot lacks bits/signum-generic.h (pulled in transitively by both).
 * Provide the few POSIX definitions we actually need instead. */
#ifndef WNOHANG
#define WNOHANG 1
#endif
#ifndef SIGTERM
#define SIGTERM 15
#endif
extern pid_t waitpid(pid_t, int *, int);
extern int kill(pid_t, int);

/* ------------------------------------------------------------------ */
/*  Internal I/O helpers                                               */
/* ------------------------------------------------------------------ */

/** Write a string to fd. Returns 0 on success, -1 on error. */
static int write_str(int fd, const char* s)
{
    size_t len = strlen(s);
    while (len > 0) {
        ssize_t n = write(fd, s, len);
        if (n < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        s   += n;
        len -= (size_t)n;
    }
    return 0;
}

/**
 * Read one line from fd into buf (up to bufsz-1 chars, NUL-terminated).
 * Strips the trailing newline.
 * Returns the line length (>= 0) or -1 on EOF/error.
 */
static int read_line(int fd, char* buf, int bufsz)
{
    int pos = 0;
    while (pos < bufsz - 1) {
        ssize_t n = read(fd, buf + pos, 1);
        if (n < 0) {
            if (errno == EINTR) continue;
            return -1;
        }
        if (n == 0) {
            /* EOF */
            if (pos == 0) return -1;
            break;
        }
        if (buf[pos] == '\n') break;
        pos++;
    }
    buf[pos] = '\0';
    /* Strip trailing \r */
    if (pos > 0 && buf[pos - 1] == '\r') {
        buf[pos - 1] = '\0';
        pos--;
    }
    return pos;
}

/**
 * Read lines from the child until [NimPulseqGUI Protocol End] is seen.
 * Accumulates the full preamble (including delimiters) into buf.
 * Returns total bytes written, or -1 on error.
 */
static int read_preamble_block(int fd, char* buf, int bufsz)
{
    int total = 0;
    char line[PULSEQLIB_BRIDGE_LINE_MAX];
    int in_block = 0;
    int len;
    int need;

    while (1) {
        len = read_line(fd, line, (int)sizeof(line));
        if (len < 0) return total > 0 ? total : -1;

        /* Append line + newline to buffer */
        need = len + 1; /* line + '\n' */
        if (total + need >= bufsz) return -1; /* buffer overflow */
        memcpy(buf + total, line, len);
        buf[total + len] = '\n';
        total += need;
        buf[total] = '\0';

        if (strstr(line, "[NimPulseqGUI Protocol]") && !in_block) {
            in_block = 1;
        }
        if (strstr(line, "[NimPulseqGUI Protocol End]")) {
            break;
        }
    }
    return total;
}

/* ------------------------------------------------------------------ */
/*  Lifecycle                                                         */
/* ------------------------------------------------------------------ */

int pulseqlib_bridge_open(pulseqlib_bridge* b,
                           const char* exe_path,
                           const char* script_path)
{
    int to_child[2];   /* host writes to_child[1], child reads to_child[0] */
    int from_child[2]; /* child writes from_child[1], host reads from_child[0] */
    pid_t pid;

    if (!b || !exe_path || !script_path) return -1;
    memset(b, 0, sizeof(*b));

    if (pipe(to_child) < 0) return -1;
    if (pipe(from_child) < 0) {
        close(to_child[0]);
        close(to_child[1]);
        return -1;
    }

    pid = fork();
    if (pid < 0) {
        close(to_child[0]);  close(to_child[1]);
        close(from_child[0]); close(from_child[1]);
        return -1;
    }

    if (pid == 0) {
        /* Child process */
        close(to_child[1]);
        close(from_child[0]);

        dup2(to_child[0], STDIN_FILENO);
        dup2(from_child[1], STDOUT_FILENO);
        close(to_child[0]);
        close(from_child[1]);

        execl(exe_path, exe_path,
              "--persistent",
              "--script", script_path,
              (char*)NULL);
        /* If we get here, exec failed */
        _exit(127);
    }

    /* Parent */
    close(to_child[0]);
    close(from_child[1]);

    b->pid        = (int)pid;
    b->to_child   = to_child[1];
    b->from_child = from_child[0];
    return 0;
}

void pulseqlib_bridge_close(pulseqlib_bridge* b)
{
    if (!b || b->pid <= 0) return;

    /* Send QUIT */
    write_str(b->to_child, "QUIT\n");

    close(b->to_child);
    close(b->from_child);
    b->to_child   = -1;
    b->from_child = -1;

    /* Wait for child with timeout: try non-blocking first */
    {
        int status;
        struct timespec ts;
        pid_t w = waitpid((pid_t)b->pid, &status, WNOHANG);
        if (w == 0) {
            /* Child still running; give it a moment then SIGTERM */
            ts.tv_sec = 0;
            ts.tv_nsec = 100000000L; /* 100 ms */
            nanosleep(&ts, NULL);
            w = waitpid((pid_t)b->pid, &status, WNOHANG);
            if (w == 0) {
                kill((pid_t)b->pid, SIGTERM);
                waitpid((pid_t)b->pid, &status, 0);
            }
        }
    }

    b->pid = 0;
}

int pulseqlib_bridge_alive(const pulseqlib_bridge* b)
{
    int status;
    if (!b || b->pid <= 0) return 0;
    if (waitpid((pid_t)b->pid, &status, WNOHANG) == 0) return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Commands                                                          */
/* ------------------------------------------------------------------ */

int pulseqlib_bridge_list_protocol(pulseqlib_bridge* b,
                                    pulseqlib_protocol* out)
{
    char preamble[4096];
    int len;
    char resp_line[PULSEQLIB_BRIDGE_LINE_MAX];

    if (!b || !out || b->pid <= 0) return -1;

    if (write_str(b->to_child, "LIST_PROTOCOL\n") < 0) return -1;

    /* Child responds: "PROTOCOL\n" then the preamble block */
    if (read_line(b->from_child, resp_line, (int)sizeof(resp_line)) < 0)
        return -1;
    if (strncmp(resp_line, "PROTOCOL", 8) != 0) return -1;

    len = read_preamble_block(b->from_child, preamble, (int)sizeof(preamble));
    if (len <= 0) return -1;

    return pulseqlib_protocol_parse(out, preamble);
}

/** Send the serialized protocol (preamble) to the child. */
static int send_protocol(pulseqlib_bridge* b,
                          const pulseqlib_protocol* proto)
{
    char buf[4096];
    int n = pulseqlib_protocol_serialize(proto, buf, (int)sizeof(buf));
    if (n <= 0) return -1;
    return write_str(b->to_child, buf);
}

int pulseqlib_bridge_validate(pulseqlib_bridge* b,
                               const pulseqlib_protocol* proto,
                               float* duration,
                               char* info, int infosz)
{
    char resp[PULSEQLIB_BRIDGE_LINE_MAX];

    if (!b || !proto || b->pid <= 0) return -1;

    if (write_str(b->to_child, "VALIDATE\n") < 0) return -1;
    if (send_protocol(b, proto) < 0) return -1;

    /* Response: "VALID <duration> <info>" or "INVALID <info>" */
    if (read_line(b->from_child, resp, (int)sizeof(resp)) < 0) return -1;

    if (strncmp(resp, "VALID ", 6) == 0) {
        /* Parse "VALID 5.32 TA = 5.32 s" */
        const char* p = resp + 6;
        if (duration) *duration = (float)atof(p);
        /* Find start of info (after the float) */
        while (*p && *p != ' ') p++;
        if (*p == ' ') p++;
        if (info && infosz > 0) {
            strncpy(info, p, infosz - 1);
            info[infosz - 1] = '\0';
        }
        return 1;
    }
    if (strncmp(resp, "INVALID", 7) == 0) {
        const char* p = resp + 7;
        if (*p == ' ') p++;
        if (info && infosz > 0) {
            strncpy(info, p, infosz - 1);
            info[infosz - 1] = '\0';
        }
        return 0;
    }

    return -1; /* unexpected response */
}

int pulseqlib_bridge_generate(pulseqlib_bridge* b,
                               const pulseqlib_protocol* proto,
                               const char* output_path)
{
    char cmd[1024];
    char resp[PULSEQLIB_BRIDGE_LINE_MAX];

    if (!b || !proto || !output_path || b->pid <= 0) return -1;

    /* "GENERATE <path>\n" */
    if (strlen(output_path) > sizeof(cmd) - 16) return -1;
    sprintf(cmd, "GENERATE %s\n", output_path);

    if (write_str(b->to_child, cmd) < 0) return -1;
    if (send_protocol(b, proto) < 0) return -1;

    /* Response: "GENERATED <path>" or "ERROR <msg>" */
    if (read_line(b->from_child, resp, (int)sizeof(resp)) < 0) return -1;

    if (strncmp(resp, "GENERATED", 9) == 0) return 0;
    return -1;
}
