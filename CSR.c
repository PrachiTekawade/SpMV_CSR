#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"

typedef struct
{
  int n_row;
  int n_col;
  int *row_ptr;
  int *col;
  double *val;
  int nnz;
}CSR;

CSR* Read_CSR(char *filename)
{
  CSR *csr = (CSR*)malloc(sizeof(CSR));
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;   
  int i, *I, *J;

  if ((f = fopen(filename, "r")) == NULL) 
    return NULL;

  if (mm_read_banner(f, &matcode) != 0)
  {
    printf("Could not process Matrix Market banner.\n");
    return NULL;
  }
  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    return NULL;

  printf("\nMatrix: %s M: %d N: %d nnz: %d\n",filename,M,N, nz);
  /* reseve memory for matrices */

  I = (int *) malloc((nz +1) * sizeof(int));
  J = (int *) malloc((nz +1) * sizeof(int));
  double* val = (double *) malloc(nz * sizeof(double));

  csr->nnz = nz;
  csr->n_row = N;
  csr->n_col = M;
  csr->row_ptr = (int*)calloc((csr->n_row )+ 2,sizeof(int));
  csr->col = (int*)calloc((csr->nnz) +1 ,sizeof(int));
  csr->val = (double*)calloc((csr->nnz) +1,sizeof(double));

  csr->row_ptr[(csr->n_row) +1] = nz+1;
  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  int pattern = mm_is_pattern(matcode);
  if(!pattern)
  {
    double temp;
    for (i=1; i<=nz; i++)
    {
      fscanf(f, "%d %d %lf\n", &I[i], &J[i],&temp);
      val[i] = temp;
    }
  }
  else
  {
    for (i=1; i<=nz; i++)
    {
      fscanf(f, "%d %d", &I[i], &J[i]);
      val[i] = 1;
    }
  }
  int index = 1;
  for (i=1;i<=csr->n_row;i++)
  {
    for (int j=1;j<=nz;j++) //iterate over I[]
    {
      if(I[j]==i)
      {
        if(csr->row_ptr[i]==0)
        {
          csr->row_ptr[i] = index;
        }
        csr->val[index] = val[j];
        csr->col[index++] = J[j];
      }
    }
  }
  fclose(f);
  free(I);
  free(J);  
  return csr;
}

int main(int argc, char *argv[]) 
{
  CSR *csr ;
  
  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    return 1;
  }
  else    
  { 
    csr = Read_CSR(argv[1]);
  }
  if(csr==NULL)
    return 1;

  double *x = (double*)calloc(((csr->n_col)+1),sizeof(double));
  double *y = (double*)calloc(((csr->n_col)+1),sizeof(double));
  for(int i=1;i<=csr->n_col;i++)
  {
    x[i]=1;
  }

  //sparse matrix vector multiplication
  clock_t start = clock();
  for(int iteration = 0 ; iteration < 1000 ; iteration++)
  {  
    for(int i=1;i<=csr->n_col;i++)
    {
      y[i] =0;
    }
    
    for(int i = 1; i<=csr->n_col;i++)
    {
      for(int j=csr->row_ptr[i];j <= csr->row_ptr[i+1]-1;j++)
      {
        y[i] = y[i] + csr->val[j] * x[csr->col[j]];
      }
    }
  }
  clock_t end = clock();

  //printing ANS to file output.txt
  FILE *out = fopen("output.txt","w");
  fprintf(out, "\nOutput vector (y):\n");
  printf("\nOutput vector (y):\n");
  for(int i=1;i<=csr->n_col;i++)
  {
    fprintf(out,"%lf\n",y[i]);
    printf("%lf\n",y[i]);
  }
  fclose(out);
  printf("\nTime required for 1000 iterations: %lf s\n", (double)(end-start)/CLOCKS_PER_SEC );
  free(y);
  free(x);
  free(csr->col);
  free(csr->row_ptr);
  free(csr->val);
  free(csr);
  return 0;
}