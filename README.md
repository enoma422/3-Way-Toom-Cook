# 3-Way-Toom-Cook

### 💻 구현 환경
|구현 환경||
|:-------------------:|:-------------------:|
| 하드웨어  | MacBook Air. Apple M2. 8GB RAM  | 
| 컴파일러 | gcc 13.1.6(-O2) |
| 정수 연산 라이브러리 | FLINT 2.9.0 |

### ⏰ 3-Way Toom-Cook 알고리듬 연산 시간 (단위: ms)
| 입력 정수의 비트 길이  | 7,680 | 15,360 |
|:-------------------:|:--------:|:--------:|
| (상황 1) | 0.006 | 0.011 |
| (상황 2) | 104.084 | 182.981 |
| (상황 3) | 184.456 | 344.708 |

### 컴파일 명령어
gcc -Wall -O2 *.c -L/opt/homebrew/Cellar/mpfr/4.1.0/lib -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -I/opt/homebrew/Cellar/mpfr/4.1.0/include -I/usr/local/include/flint/ -lflint -lmpfr -lgmp -lpthread
