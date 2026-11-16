// prof_asm.c -- SIGPROF sampling profiler for aarch64 static binaries.
// Records PCs into a fixed buffer; dumps them on request. Enabled only when
// the caller decides (env var PROFILE_US).
#define _GNU_SOURCE
#include <signal.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ucontext.h>

#define PROFSZ (1 << 21)
static volatile uint64_t prof_pcs[PROFSZ];
static volatile uint64_t prof_n;

static void prof_handler(int sig, siginfo_t *si, void *uc)
{
    (void)sig; (void)si;
    ucontext_t *u = (ucontext_t *)uc;
    uint64_t pc = (uint64_t)u->uc_mcontext.pc;
    uint64_t n = prof_n;
    if (n < PROFSZ) {
        prof_pcs[n] = pc;
        prof_n = n + 1;
    }
}

void prof_start(int usec)
{
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_sigaction = prof_handler;
    sa.sa_flags = SA_SIGINFO | SA_RESTART;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGPROF, &sa, NULL);
    struct itimerval it;
    it.it_interval.tv_sec = 0;
    it.it_interval.tv_usec = usec;
    it.it_value = it.it_interval;
    setitimer(ITIMER_PROF, &it, NULL);
}

void prof_stop(void)
{
    struct itimerval it;
    memset(&it, 0, sizeof it);
    setitimer(ITIMER_PROF, &it, NULL);
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
