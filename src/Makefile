CFLAGS=-g -c -Wall -Wextra -Werror -std=c11 -pedantic
FLAGS=-Wall -Wextra -Werror -std=c11 -pedantic
ASAN=#-fsanitize=address -fanalyzer


all: s21_matrix.a

s21_matrix.a:
	gcc $(CFLAGS) -o s21_matrix.o -c s21_matrix.c
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

s21_matrix.a_for_tests:
	rm -f s21_s21_matrix.a
	gcc $(CFLAGS) $(ASAN) -fprofile-arcs -ftest-coverage -o s21_matrix.o -c s21_matrix.c
	ar rc s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

clean:
	rm -f *.a *.o
	rm -f *.info *.gcda *.gcno *.gcov
	rm -rf coverage
	rm -f mytests

test: s21_matrix.a_for_tests
	gcc $(FLAGS) $(ASAN) -fprofile-arcs -ftest-coverage ./tests/s21_test.c -o mytests -lcheck -L. s21_matrix.a
	./mytests
	geninfo . -b . -o ./coverage.info
	genhtml coverage.info -o ./coverage

gcov_report:
	open coverage/index.html
