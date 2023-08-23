#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int code = OK;

  if (rows < 1 || columns < 1) {
    code = MATRIX_ERROR;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    if (result->matrix) {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = calloc(columns, sizeof(double));
      }
    } else {
      code = MATRIX_ERROR;
    }
  }

  return code;
}

int s21_matrix_exist(matrix_t *A) {
  int code = OK;

  if (!A->matrix || A->columns < 1 || A->rows < 1) {
    code = MATRIX_ERROR;
  }

  return code;
}

void s21_remove_matrix(matrix_t *A) {
  if (s21_matrix_exist(A) == OK) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int code = SUCCESS;

  if (s21_matrix_exist(A) == OK && s21_matrix_exist(B) == OK &&
      A->rows == B->rows && A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if ((A->matrix[i][j] > 0 && B->matrix[i][j] < 0) ||
            (A->matrix[i][j] < 0 && B->matrix[i][j] > 0) ||
            (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7)) {
          code = FAILURE;
        }
      }
    }
  } else {
    code = FAILURE;
  }

  return code;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK && s21_matrix_exist(B) == OK) {
    if (A->rows == B->rows && A->columns == B->columns &&
        s21_create_matrix(A->rows, B->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK && s21_matrix_exist(B) == OK) {
    if (A->rows == B->rows && A->columns == B->columns &&
        s21_create_matrix(A->rows, B->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK) {
    if (!isnan(number) &&
        s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK && s21_matrix_exist(B) == OK) {
    if (A->columns == B->rows &&
        s21_create_matrix(A->rows, B->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int k = 0; k < A->columns; k++) {
          for (int j = 0; j < B->columns; j++) {
            result->matrix[i][j] =
                result->matrix[i][j] + A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK) {
    if (s21_create_matrix(A->columns, A->rows, result) == OK) {
      for (int i = 0; i < A->columns; i++) {
        for (int j = 0; j < A->rows; j++) {
          result->matrix[i][j] = A->matrix[j][i];
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_minor_matrix(matrix_t *A, int target_row, int target_column,
                     matrix_t *minor_result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK) {
    if (A->rows > 1 &&
        s21_create_matrix(A->rows - 1, A->columns - 1, minor_result) == OK) {
      int k = 0, m = 0;
      for (int i = 0; i < A->rows; i++) {
        if (i == target_row) continue;
        for (int j = 0; j < A->columns; j++) {
          if (j == target_column) continue;
          minor_result->matrix[k][m++] = A->matrix[i][j];
        }
        k++;
        m = 0;
      }
    } else if (A->rows == 1 && A->columns == 1 &&
               s21_create_matrix(A->rows, A->columns, minor_result) == OK) {
      minor_result->matrix[0][0] = A->matrix[0][0];
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

void s21_determinant2x2(matrix_t *A, double *result) {
  if (A->rows == 2 && A->columns == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  }
}

int s21_determinant(matrix_t *A, double *result) {
  int code = OK;
  *result = 0;

  if (s21_matrix_exist(A) == OK) {
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        s21_determinant2x2(A, result);
      } else {
        for (int j = 0; j < A->rows; j++) {
          double tmp_det = 0;
          matrix_t tmp = {0};
          s21_minor_matrix(A, 0, j, &tmp);
          s21_determinant(&tmp, &tmp_det);
          *result += pow(-1, j) * (A->matrix[0][j] * tmp_det);
          s21_remove_matrix(&tmp);
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int code = OK;

  if (s21_matrix_exist(A) == OK) {
    if (A->rows == A->columns &&
        s21_create_matrix(A->columns, A->rows, result) == OK) {
      for (int i = 0; i < A->rows && code == OK; i++) {
        for (int j = 0; j < A->columns && code == OK; j++) {
          double det = 0;
          matrix_t tmp = {0};
          if (s21_minor_matrix(A, i, j, &tmp) == OK &&
              s21_determinant(&tmp, &det) == OK) {
            result->matrix[i][j] = pow(-1, i + j) * det;
            s21_remove_matrix(&tmp);
          } else {
            code = CALCULATION_ERROR;
          }
        }
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int code = OK;
  double det = 0;

  if (s21_matrix_exist(A) == OK) {
    if (s21_determinant(A, &det) == OK && det != 0) {
      if (A->rows == 1 &&
          s21_create_matrix(A->rows, A->columns, result) == OK) {
        result->matrix[0][0] = 1 / A->matrix[0][0];
      } else {
        matrix_t tmp = {0}, mult_tmp = {0};
        s21_calc_complements(A, &tmp);
        s21_transpose(&tmp, &mult_tmp);
        s21_mult_number(&mult_tmp, 1 / det, result);
        s21_remove_matrix(&tmp);
        s21_remove_matrix(&mult_tmp);
      }
    } else {
      code = CALCULATION_ERROR;
    }
  } else {
    code = MATRIX_ERROR;
  }

  return code;
}

void s21_fill_matrix(double num, matrix_t *result) {
  if (s21_matrix_exist(result) == OK) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = num++;
      }
    }
  }
}
