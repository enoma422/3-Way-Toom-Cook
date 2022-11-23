// ---------------------------------------------
// Toom-Cook 3 Way multiplication
// ---------------------------------------------
void mul_TC(fmpz_t out, fmpz_t a, fmpz_t b){
    slong n = (fmpz_bits(a) / w) + 1;
    slong m = (fmpz_bits(b) / w) + 1;
    
    if(10 >= MIN(n, m)){
        fmpz_mul(out, a, b);
        return;
    }

    slong l = (MAX(n, m) + 1) / 3;
    slong lw = l * w;
    slong lw2 = lw * 2;
    
    // --------------------------------------------------------------------------------------------
    // Step 1: splitting
    // --------------------------------------------------------------------------------------------
    // input a, b를 3파트로 분할
    // a = a2 * 2^{lw * 2} + a1 * 2^{lw} + a0
    // b = b2 * 2^{lw * 2} + b1 * 2^{lw} + b0
    fmpz_t a2, a1, a0;
    fmpz_t b2, b1, b0;
    fmpz_init(a2);    fmpz_init(a1);    fmpz_init(a0);
    fmpz_init(b2);    fmpz_init(b1);    fmpz_init(b0);

    fmpz_fdiv_q_2exp(a2, a, lw2);    fmpz_fdiv_q_2exp(a1, a, lw);
    fmpz_fdiv_r_2exp(a1, a1, lw);    fmpz_fdiv_r_2exp(a0, a, lw);
    fmpz_fdiv_q_2exp(b2, b, lw2);    fmpz_fdiv_q_2exp(b1, b, lw);
    fmpz_fdiv_r_2exp(b1, b1, lw);    fmpz_fdiv_r_2exp(b0, b, lw);

    // --------------------------------------------------------------------------------------------
    // Step 2: Evaluation
    // --------------------------------------------------------------------------------------------
    fmpz_t p0, p1, pm1, pm2, pinf;
    fmpz_t q0, q1, qm1, qm2, qinf;
    fmpz_init(p0);    fmpz_init(p1);    fmpz_init(pm1);    fmpz_init(pm2);    fmpz_init(pinf);
    fmpz_init(q0);    fmpz_init(q1);    fmpz_init(qm1);    fmpz_init(qm2);    fmpz_init(qinf);

    // p0 = a0
    fmpz_init_set(p0, a0);
    // p1 = a0 + a1 + a2
    // pm1 = a0 - a1 + a2
    fmpz_add(p1, a0, a2);    fmpz_sub(pm1, p1, a1);    fmpz_add(p1, p1, a1);
    // pm2 = a0 - 2a1 + 4a2
    fmpz_init_set(pm2, a0);    fmpz_addmul_ui(pm2, a2, 4);    fmpz_submul_ui(pm2, a1, 2);
    // pinf = a2
    fmpz_init_set(pinf, (a2));

    fmpz_init_set(q0, (b0));
    fmpz_add(q1, b0, b2);    fmpz_sub(qm1, q1, b1);    fmpz_add(q1, q1, b1);
    fmpz_init_set(qm2, (b0));    fmpz_addmul_ui(qm2, b2, 4);    fmpz_submul_ui(qm2, b1, 2);
    fmpz_init_set(qinf, (b2));

    // --------------------------------------------------------------------------------------------
    // Step 3: Pointwise multiplication
    // --------------------------------------------------------------------------------------------
    fmpz_t r0, r1, rm1, rm2, rinf;
    fmpz_init(r0);    fmpz_init(r1);    fmpz_init(rm1);    fmpz_init(rm2);    fmpz_init(rinf);
    mul_TC(r0, p0, q0);
    mul_TC(r1, p1, q1);
    mul_TC(rm1, pm1, qm1);
    mul_TC(rm2, pm2, qm2);
    mul_TC(rinf, pinf, qinf);
    
    // --------------------------------------------------------------------------------------------
    // Step 4: Interpolation
    // --------------------------------------------------------------------------------------------
    fmpz_t rx0, rx1, rx2, rx3, rx4;
    fmpz_init(rx0);    fmpz_init(rx1);    fmpz_init(rx2);    fmpz_init(rx3);    fmpz_init(rx4);
    
    fmpz_init_set(rx0, r0);
    
    fmpz_init_set(rx4, rinf);
    
    fmpz_sub(rx3, rm2, r1);
    fmpz_cdiv_q_si(rx3, rx3, 3);
    
    fmpz_sub(rx1, r1, rm1);
    fmpz_cdiv_q_si(rx1, rx1, 2);
    
    fmpz_sub(rx2, rm1, r0);

    fmpz_sub(rx3, rx2, rx3);
    fmpz_cdiv_q_si(rx3, rx3, 2);   
    fmpz_addmul_si(rx3, rinf, 2);

    fmpz_add(rx2, rx2, rx1);
    fmpz_sub(rx2, rx2, rx4);

    fmpz_sub(rx1, rx1, rx3);    
    
    // --------------------------------------------------------------------------------------------
    // Step 5: recomposition
    // --------------------------------------------------------------------------------------------
    fmpz_mul_2exp(rx1, rx1, lw);
    fmpz_mul_2exp(rx2, rx2, lw*2);
    fmpz_mul_2exp(rx3, rx3, lw*3);
    fmpz_mul_2exp(rx4, rx4, lw*4);

    fmpz_add(out, out, rx0);
    fmpz_add(out, out, rx1);
    fmpz_add(out, out, rx2);
    fmpz_add(out, out, rx3);
    fmpz_add(out, out, rx4);

    // --------------------------------------------------------------------------------------------
    fmpz_clear(a2);    fmpz_clear(a1);    fmpz_clear(a0);
    fmpz_clear(b2);    fmpz_clear(b1);    fmpz_clear(b0);
    fmpz_clear(p0);    fmpz_clear(p1);    fmpz_clear(pm1);    fmpz_clear(pm2);    fmpz_clear(pinf);
    fmpz_clear(q0);    fmpz_clear(q1);    fmpz_clear(qm1);    fmpz_clear(qm2);    fmpz_clear(qinf);
    fmpz_clear(r0);    fmpz_clear(r1);    fmpz_clear(rm1);    fmpz_clear(rm2);    fmpz_clear(rinf);

}


int check_mul(fmpz_t a, fmpz_t b, fmpz_t res){
    fmpz_t tmp;
    fmpz_init(tmp);
    fmpz_mul(tmp, a, b);
  
    // --------------------------------------------------------
    // Equal
    // --------------------------------------------------------
    if(fmpz_equal(tmp, res)){
        fmpz_clear(tmp);
        // printf(" equal");
        return 0;
    }
    // --------------------------------------------------------
    // Not Equal
    // --------------------------------------------------------
    printf("\n a : ");
    fmpz_print(a);
    printf("\n b : ");
    fmpz_print(b);
    printf("\n----- 일반 곱셈 -----\n");
    printf(" a * b : ");
    fmpz_print(tmp);
    printf("\n----- 3-Way Toom-Cook 곱셈 -----\n");
    printf(" a * b : ");
    fmpz_print(res);

    fmpz_clear(tmp);
    return 1;
}


int main(){
    // -------------------
    // 반복 횟수
    // -------------------
    // int cnt = 1;
    // int cnt = 15360;
    int cnt = 30720;
    
    // -------------------
    // 입력되는 두 정수의 비트 길이
    // -------------------
    int bitlen = 7680;
    // int bitlen = 15360;
    
    fmpz_t test_a;    fmpz_t test_b;
    fmpz_init(test_a);      fmpz_init(test_b);
        
    FLINT_TEST_INIT(state);
    double time_res = 0;

    for(int i=0; i<cnt; i++){    
        fmpz_t test_res;
        fmpz_init(test_res);
        fmpz_randtest(test_a, state, bitlen);
        fmpz_randtest(test_b, state, bitlen);

        clock_t start = clock();
  
        // 3-Way Toom-Cook Multiplication
        mul_TC(test_res, test_a, test_b);

        clock_t end = clock();

        time_res += (double) (end - start);
        if(check_mul(test_a, test_b, test_res)) {
            printf("\nerror");
        }
        fmpz_clear(test_res);
    }
    printf(" %d회 반복 >> %lf 밀리초\n", cnt, (time_res * 1000) / (CLOCKS_PER_SEC));
  
    fmpz_clear(test_a);
    fmpz_clear(test_b);
    
    return 0;
}
