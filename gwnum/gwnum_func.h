#pragma once

typedef int (*c_func1r)(struct gwasm_data *);
typedef void (*c_func1_uptr)(struct gwasm_data *, uintptr_t);

void pass1_wake_up_threads(struct gwasm_data *asm_data);
void pass1_pre_carries(struct gwasm_data *asm_data);
int pass1_post_carries(struct gwasm_data *asm_data);
//int pass1_get_next_block_mt(struct gwasm_data *asm_data);
void pass2_wake_up_threads(struct gwasm_data *asm_data);
//int pass2_get_next_block_mt(struct gwasm_data *asm_data);

void zpad_sub7(struct gwasm_data *asm_data);
