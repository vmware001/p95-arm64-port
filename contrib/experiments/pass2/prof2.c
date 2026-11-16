// prof2.c -- high-resolution SIGPROF PC sampling profiler for aarch64 static binaries.
// Uses POSIX interval timer on CLOCK_MONOTONIC (hrtimer, ~us resolution) instead of
// setitimer/ITIMER_PROF which is limited to 1/HZ jiffies (4 ms on this kernel).
#define _GNU_SOURCE
#include <signal.h>
#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ucontext.h>

#define PROFSZ (1 << 21)
static volatile uint64_t prof_pcs[PROFSZ];
static volatile uint64_t prof_n;
static timer_t prof_tid;

static void prof_handler(int sig, siginfo_t *si, void *uc)
{
    (void)sig; (void)si;
    ucontext_t *u = (ucontext_t *)uc;
    uint64_t pc = (uint64_t)u->uc_mcontext.pc;
    uint64_t n = prof_n;
    if (n < PROFSZ) { prof_pcs[n] = pc; prof_n = n + 1; }
}

void prof_start(int usec)
{
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_sigaction = prof_handler;
    sa.sa_flags = SA_SIGINFO | SA_RESTART;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGPROF, &sa, NULL);
    struct sigevent sev;
    memset(&sev, 0, sizeof sev);
    sev.sigev_notify = SIGEV_SIGNAL;
    sev.sigev_signo = SIGPROF;
    timer_create(CLOCK_MONOTONIC, &sev, &prof_tid);
    struct itimerspec its;
    memset(&its, 0, sizeof its);
    its.it_interval.tv_sec = usec / 1000000;
    its.it_interval.tv_nsec = (long)(usec % 1000000) * 1000L;
    its.it_value = its.it_interval;
    timer_settime(prof_tid, 0, &its, NULL);
}

void prof_stop(void)
{
    struct itimerspec z;
    memset(&z, 0, sizeof z);
    timer_settime(prof_tid, 0, &z, NULL);
    timer_delete(prof_tid);
    signal(SIGPROF, SIG_IGN);
}

void prof_dump(const char *label)
{
    uint64_t n = prof_n;
    fprintf(stderr, "PROF %s n=%llu\n", label, (unsigned long long)n);
    for (uint64_t i = 0; i < n; i++)
        fprintf(stderr, "PC %llx\n", (unsigned long long)prof_pcs[i]);
    fflush(stderr);
}

void prof_reset(void) { prof_n = 0; }
