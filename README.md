# 3-Way-Toom-Cook

### 💻 구현 환경

| 하드웨어  | MacBook Air. Apple M2. 8GB RAM  | 
|:-------------------:|:-------------------:|
| **컴파일러** | **gcc 13.1.6(-O2)** |
| **정수 연산 라이브러리** | **FLINT 2.9.0** |

### ⏰ 3-Way Toom-Cook 알고리듬 연산 시간 (단위: ms)
| 입력 정수의 비트 길이  | 7,680 | 15,360 |
|:-------------------:|:--------:|:--------:|
| (상황 1) | 0.006 | 0.011 |
| (상황 2) | 104.084 | 182.981 |
| (상황 3) | 184.456 | 344.708 |

### 컴파일 명령어
gcc -Wall -O2 *.c -L/opt/homebrew/Cellar/mpfr/4.1.0/lib -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -I/opt/homebrew/Cellar/mpfr/4.1.0/include -I/usr/local/include/flint/ -lflint -lmpfr -lgmp -lpthread*

### 3-Way-Toom-Cook 유사 부호
```
procedure MUL^{3ToomCook}\left(A,B\right)
If flag\geq min\left(len\left(A\right),len\left(B\right)\right) then
	return AB	\vartriangleright AB:일반곱셈
end if
l\gets max\left(len\left(A\right),len\left(B\right)\right)\ /\ 3 
A_2,{\ A}_1,\ {\ A}_0\gets A_{\left[n:2l\right]},\ \ A_{\left[2l:l\right]},\ \ A_{\left[l:0\right]} 
B_2,{\ B}_1,\ {\ B}_0\gets B_{\left[n:2l\right]},\ \ B_{\left[2l:l\right]},\ \ B_{\left[l:0\right]} 
cnt\gets0 
T\gets[\ ] 
for i in \left[-2,-1,0,1,\infty\right] do
	T[cnt]←MUL3ToomCookpi,qi 

	cnt\gets cnt+1 
r\gets T\cdot\mathcal{M}^T 
R_{1\left[6l:4l\right]},R_{1\left[4l:2l\right]},R_{1\left[2l:0\right]}\gets r\left[4\right],r\left[2\right],r\left[0\right]  
R_{0\left[5l:3l\right]},R_{0\left[3l:l\right]},R_{0\left[l:0\right]}\gets r\left[3\right],r\left[1\right],0_l 
R\gets R_1+R_0 
return R
end procedure

```
