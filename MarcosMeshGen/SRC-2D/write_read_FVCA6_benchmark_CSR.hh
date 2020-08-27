/**************************************************************************
 
  The solver-comp module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-common and dune-istl modules 
  providing a possibility to compare run times of several linear
  solvers for different discretization schemes. 

  Copyright (C) 2010 Robert Kloefkorn

  The benchmark-runtime module is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.

**************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* DO NOT CHANGE THIS VARIABLE, YOU'LL GET WHAT YOU DESERVE */
static const double SOLVERBENCH_VERSION = 0.67 ; 
static const double ACCEPTED_VERSION = 0.58 ;

#ifndef WRITEXDR_H_INCLUDED
#define WRITEXDR_H_INCLUDED

/* system headers */
#include <assert.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

/**
     @defgroup writeXDR Write Benchmark files using XDR
 */

/**
     @defgroup readXDR Reading the Solution using XDR
 */

typedef void rw_int_t(int *);
typedef void rw_double_t(double *);
typedef void rw_string_t(char *, const unsigned int );

/* XDR write and read routines */
void readwrite_int_xdr(int * value);
void readwrite_double_xdr(double * value);
void readwrite_string_xdr(char* value, const unsigned int );

/** \ingroup writeXDR 
    \author Robert Kloefkorn
    \brief Implementation of the write routine for calls from \b C/C++ \b programs.
    The implementation writes the number of unknowns, the matrix, and the right hand
    side which then can be read with the program `benchruntime' 

    \param[in]  path             path to write file to, e.g. "./"
    \param[in]  nameAndScheme    contributor and name of the scheme 
    \param[in]  test             number of the test, valid are { 1,..., 5 }
    \param[in]  mesh             name of the mesh file used for calculation, e.g. "tet.2.msh" or "vmesh_4.msh"
    \param[in]  format           format of CSR storage, values are {0,...,3} \n
                                  \b 0 = offset in rows and columns is zero (default \b C \b format) \n
                                  \b 1 = offset in rows is 0, in columns 1 \n 
                                  \b 2 = offset in rows is 1, in columns 0 \n 
                                  \b 3 = offset in rows and columns is 1 (default \b Fortran \b format)
    \param[in]  blockSize        block size of DoFs ( \b = \b 1 for most schems such as Finite Volume and Finite Element schemes) 
    \param[in]  numberOfRows     number of rows of the matrix 
    \param[in]  numberOfColumns  number of columns of the matrix  
    \param[in]  rows             offset for each row of the matrix  [0,...,numberOfRows], note length is numberOfRows+1, 
                                 the last entries is the number of non-zeros + offset
    \param[in]  columns          real columns number of each matrix entry  [0,...,numberNonZeros-1], length is numberNonZeros 
    \param[in]  matrix           matrix entiries [0,...,numberNonZeros-1], length is numberNonZeros
    \param[in]  rhs              right hand side [0,...,numberOfRows-1], length is numberOfRows

    \note Only call this routine from \b C/C++ \b programs.\n
          An example can be found in \ref example_C.c.
    
*/
void writeBenchmarkFile(const char* path,           /* path to write file to, e.g. `./' */
                        const char* nameAndScheme,  /* contributor and scheme */
                        const int test,             /* number of test , possible set { 1,..., 5 } */
                        const char* mesh,           /* mesh file used for calculation */
                        const int format,           /* format of CSR storage */
                        const int blockSize,        /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                        const int numberOfRows,     /* number of rows of the matrix  */
                        const int numberOfColumns,  /* number of columns of the matrix */
                        const int* rows,            /* the offset for each row in the vector matrix */
                        const int* columns,         /* the columns number for each matrix entry */
                        const double* matrix,       /* the matrix entries */
                        const double* rhs);         /* the right hand side */ 

/** \ingroup writeXDR 
    \author Robert Kloefkorn
    \brief Implementation of the write routine for calls from \b Fortran \b programs.
    The implementation writes the number of unknowns, the matrix, and the right hand
    side which then can be read with the program `benchruntime' 

    \param[in]  path             path to write file to, e.g. "./"
    \param[in]  nameAndScheme    contributor and name of the scheme 
    \param[in]  test             number of the test, valid are { 1,..., 5 }
    \param[in]  mesh             name of the mesh file used for calculation, e.g. "tet.2.msh" or "vmesh_4.msh"
    \param[in]  format           format of CSR storage, values are [0,...,3] \n
                                  \b 0 = offset in rows and columns is zero (default \b C \b format) \n
                                  \b 1 = offset in rows is 0, in columns 1 \n 
                                  \b 2 = offset in rows is 1, in columns 0 \n 
                                  \b 3 = offset in rows and columns is 1 (default \b Fortran \b format)
    \param[in]  blockSize        block size of DoFs ( \b = \b 1 for most schems such as Finite Volume and Finite Element schemes) 
    \param[in]  numberOfRows     number of rows of the matrix 
    \param[in]  numberOfColumns  number of columns of the matrix  
    \param[in]  rows             offset for each row of the matrix  [0,...,numberOfRows], note length is numberOfRows+1, the last entries is the number of non-zeros + offset
    \param[in]  columns          real columns number of each matrix entry  [0,...,numberNonZeros-1], length is numberNonZeros 
    \param[in]  matrix           matrix entiries [0,...,numberNonZeros-1], length is numberNonZeros
    \param[in]  rhs              right hand side [0,...,numberOfRows-1], length is numberOfRows

    \note Only call this routine from \b Fortran \b programs. \n
          The call look as follows:  \n
          \b call \b writebenchmarkfile( CSRformat, filename, blocksize, n, n, rows, columns, matrix, rhs ) \n
          An example can be found in \ref example_F.f90.
*/
void writebenchmarkfile_(const char* path,           /* path to write file to */
                         const char* nameAndScheme,  /* filename of file to write to */
                         const int* test,            /* number of test , possible set { 1,..., 5 } */
                         const char* mesh,           /* mesh file used for calculation */
                         const int* format,          /* format of CSR storage */
                         const int* blockSize,       /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                         const int* numberOfRows,    /* number of rows of the matrix  */
                         const int* numberOfColumns, /* number of columns of the matrix */
                         const int* rows,            /* the offset for each row in the vector matrix */
                         const int* columns,         /* the columns number for each matrix entry */
                         const double* matrix,       /* the matrix entries */
                         const double* rhs);         /* the right hand side */ 


/** \ingroup readXDR 
    \author Robert Kloefkorn
    \brief Implementation of the read routine 
      for calls from \b C/C++ \b programs 

    \param[in]  filename          file name that contains the solution 
    \param[in]  numberOfUknowns   number of unknowns to be read (for checking) 
    \param[out] solution          pointer to memory where the solution is stored which already has the
                                  correct length given by the discretization grid 
           
    \note Only call this routine from \b C/C++ \b programs.\n
          An example can be found in \ref read_C.c.
    
 */
void readSolution(const char* filename,  /* filename to read solution from */
                  const int numberOfUknowns, /* number of unknowns to be read (for checking) */
                  double* solution);  /* pointer to memory for solution */

/** \ingroup readXDR 
    \author Robert Kloefkorn
    \brief Implementation of the read routine 
           for calls from \b Fortran \b programs 

    \param[in]  filename          file name that contains the solution 
    \param[in]  numberOfUknowns   number of unknowns to be read (for checking) 
    \param[out] solution          pointer to memory where the solution is stored which already has the
                                  correct length given by the discretization grid 

    \note Only call this routine from \b Fortran \b programs. \n
          The call look as follows:  \n
          \b call \b readsolution( filename, n, solution ) \n
          An example can be found in \ref read_F.f90.
 */
void readsolution_(const char* filename,  /* filename to read solution from */
                   const int* numberOfUknowns, /* number of unknowns to be read (for checking) */
                   double* solution); /* the solution */

/***************************************************************************
 *
 *  Helper functions 
 *
 **************************************************************************/

/* read-write of dimensions */
int  readWriteHeader(rw_int_t* rw_int,      /* read-write integer function pointer */
                     rw_double_t* rw_double,/* read-write double function pointer */
                     rw_string_t* rw_string,/* read-write strings function pointer */
                     char** nameAndScheme,  /* filename of file to write to */
                     int* test,             /* number of test , possible set { 1,..., 5 } */
                     char** mesh,           /* mesh file name */
                     int* blockSize,        /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                     int* numberNonZeros,   /* number of non-zero entries of the matrix */
                     int* numberOfRows,     /* number of rows of the matrix  */
                     int* numberOfColumns); /* number of columns of the matrix */

/* read-write matrix enrites */
void 
readWriteMatrixAndRhs(rw_int_t* rw_int,    /*  read-write integer function pointer */
                      rw_double_t* rw_double, /* read-write double function pointer */
                      const int numberOfRows,    /* number of rows of the matrix  */
                      const int numberOfColumns, /* number of columns of the matrix */
                      int* rows,                 /* the offset for each row in the vector matrix */
                      int* columns,              /* the columns number for each matrix entry */
                      double* matrix,            /* the matrix entries */
                      double* rhs,               /* the right hand side */ 
                      const int coloffset);      /* offset for indices according to format */

/* read-write solution  */
int readWriteSolution(rw_int_t* rw_int,
                      rw_double_t* rw_double,
                      const int numberOfUnknowns,  /* number of rows of the matrix  */
                      double* solution);           /* the solution */

/* read only block size from given file */
int readBlockSize(const char* filename);   /* filename to read blocksize from */

/* set global xdr pointer (write mode) */
void setXDRWritePointer( XDR* xdrs );
/* set global xdr pointer (read mode) */
void setXDRReadPointer( XDR* xdrs );
#endif



XDR* xdrPtr = NULL;
int writeMode = 0;

/* offset for rows/columns
 *
 * Format is 0 = (0,0) 
 *           1 = (0,1) 
 *           2 = (1,0) 
 *           3 = (1,1) 
 */
int getColumnOffSet( const int format, 
                     const int numberOfRows,
                     const int* rows, 
                     const int* columns)
{
  const int offsets[4][2] = { {0,0}, {0,1}, {1,0}, {1,1} };
  const int rOffset = rows[ 0 ];
  const int nnz = rows[ numberOfRows ] - rOffset;
  int i;
  int cOffset = -1;

  if( format < 0 || format > 3 ) 
  {
    fprintf(stderr,"ERROR: wrong format (%d) given, should be out of 0,...,3 !\n",format);
    exit(EXIT_FAILURE);
  }

  if( offsets[ format ][ 0 ] != rOffset )
  {
    fprintf(stderr,"ERROR: wrong row offset (%d) given, should be %d !\n",rOffset, offsets[ format ][ 0 ] );
    exit(EXIT_FAILURE); 
  }
  
  for( i=0; i<nnz; ++i )
  {
    if( columns[ i ] == offsets[ format ][ 1 ] )
    {
      cOffset = offsets[ format ][ 1 ];
      break ;
    }
  }

  if( offsets[ format ][ 1 ] != cOffset )
  {
    fprintf(stderr,"WARNING: couldn't determine columns offset, using %d instead !", offsets[ format ][ 1 ] );
    return offsets[ format ][ 1 ] ;
  }

  return cOffset;
}


/* XDR write and read routines */
void readwrite_int_xdr(int * value)
{
  assert( xdrPtr );
  xdr_int( xdrPtr, value);
}

/* XDR write and read routines */
void readwrite_double_xdr(double * value)
{
  assert( xdrPtr );
  xdr_double( xdrPtr, value);
}

/* XDR write and read routines for strings */
void readwrite_string_xdr(char * value, const unsigned int len )
{
  assert( xdrPtr );
  xdr_string( xdrPtr, &value, len );
}

/* read and write version of package */
void readWriteVersion( rw_double_t* rw_double) 
{
  double version = SOLVERBENCH_VERSION;
  /* check version */
  rw_double( &version );
  if( version < ACCEPTED_VERSION ) 
  {
    fprintf(stderr,"ERROR: Version of file is (%lf), we need at least (%lf)! \n",version,ACCEPTED_VERSION);
    exit( EXIT_FAILURE );
  }
}

/* read-write of dimensions */
int  readWriteHeader(rw_int_t* rw,          /* read-write function pointer */
                     rw_double_t* rw_double,/* read-write double function pointer */
                     rw_string_t* rw_string,/* read-write function point for strings */
                     char** nameAndScheme,  /* filename of file to write to */
                     int* test,             /* number of test , possible set { 1,..., 5 } */
                     char** mesh,           /* mesh file name */
                     int* blockSize,        /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                     int* numberNonZeros,   /* number of non-zero entries of the matrix */
                     int* numberOfRows,     /* number of rows of the matrix  */
                     int* numberOfColumns)  /* number of columns of the matrix */
{
  int check ;
  const char* fvca6 = "BENCHMARKFVCA6";
  const int length = strlen( fvca6 ) + 16;
  char* fvcaStr = (char *) malloc( length * sizeof( char ));
  int len = 0;
  int freeMem = 0;

  /* check version of file */
  readWriteVersion( rw_double );

  assert( fvcaStr );
  /* set string only  in write mode */
  if ( writeMode )
  {
    sprintf( fvcaStr,"%s", fvca6 );
  }
  else 
  {
    /* in read mode set some string to have correct length */
    sprintf( fvcaStr,"FVCA6BENCHMARK");
  }

  /* store or read dimensions */
  rw(blockSize);

  rw(numberNonZeros);
  rw(numberOfRows);
  rw(numberOfColumns);

  rw( test );

  rw_string( fvcaStr, strlen( fvcaStr) );

  if( writeMode )
    len = strlen( *nameAndScheme );

  rw( &len );
  if( *nameAndScheme == NULL ) 
  {
    *nameAndScheme = (char *) malloc( (len + 16) * sizeof( char ) ); 
    freeMem = 1;
  }

  rw_string( *nameAndScheme, len );

  if( writeMode )
    len = strlen( *mesh );

  rw( &len );
  if( *mesh == NULL ) 
  {
    assert( freeMem );
    *mesh = (char *) malloc( (len+16) * sizeof( char ) ); 
  }
  rw_string( *mesh , len );

  check = strcmp( fvca6, fvcaStr );
  free( fvcaStr );

  if( check != 0 ) 
  {
    fprintf(stderr,"ERROR: string `%s' not found in xdr file!\n", fvca6);
    exit(1);
  }
  return freeMem;
}

/* read-write matrix enrites */
void readWriteMatrixAndRhs(
                          rw_int_t* rw_int,
                          rw_double_t* rw_double,
                          const int numberOfRows,    /* number of rows of the matrix  */
                          const int numberOfColumns, /* number of columns of the matrix */
                          int* rows,            /* the offset for each row in the vector matrix */
                          int* columns,         /* the columns number for each matrix entry */
                          double* matrix,       /* the matrix entries */
                          double* rhs,          /* the right hand side */ 
                          const int offset)     /* offset for indices (0 for C and C++, 1 for Fortran) */
{
  int i, val;
  const int rowoffset = rows[ 0 ]; 

  int numberNonZeros = 0;

  /* read-write all offsets of rows */
  for(i=0; i<numberOfRows + 1; ++i) 
  {
    /* apply possible row offset */
    val = rows[ i ] - rowoffset;
    rw_int( &val );
    /* re-adjust rows numbering */
    rows[ i ] = val + rowoffset;

#ifdef PRINT_OUTPUT
    printf("%d \n", val);
#endif

  }

  /* set number of non zeros */
  /* NOTE: this has to be done here because on reading rows have to be
   * read first */
  numberNonZeros = rows[ numberOfRows ] - rowoffset;

#ifdef PRINT_OUTPUT
  printf("%d %d %d \n", numberOfRows, numberOfColumns, numberNonZeros );
  printf("*********************************\n");
#endif

  /* read-write all entries and columns numbers */
  for(i=0; i<numberNonZeros; ++i) 
  {
    /* apply offset */
    val = columns[i] - offset ;
    rw_int( &val );
    /* for reading */
    columns[ i ] = val + offset;

#ifdef PRINT_OUTPUT
    printf("%d   %d \n", val,columns[i]);
#endif

  }

#ifdef PRINT_OUTPUT
  printf("*********************************\n");
#endif

  /* read-write all entries and columns numbers */
  for(i=0; i<numberNonZeros; ++i) 
  {
    rw_double( &matrix[i] );
    /*
    if( fabs( matrix[i] ) < 1e-14 )
    {
      fprintf(stderr,"ERROR: matrix entry (%d) is very small %f (almost zero), store only non-zero entries.\n", i, matrix[ i ]);
    }
    */
#ifdef PRINT_OUTPUT
    printf("%f \n", matrix[i]);
#endif
  }

#ifdef PRINT_OUTPUT
  printf("*********************************\n");
#endif

  /* read-write right hand side */
  for(i=0; i<numberOfRows; ++i) 
  {
    rw_double( &rhs[i] );
#ifdef PRINT_OUTPUT
    printf("%f \n",rhs[i]);
#endif
  }

#ifdef PRINT_OUTPUT
  printf("*********************************\n");
#endif
}


void doWriteBenchmarkFile(const char* filename,      /* filename of file to write to */
                          const char* nameAndScheme, /* filename of file to write to */
                          const int test,            /* number of test , possible set { 1,..., 5 } */
                          const char* mesh,          /* mesh file name */ 
                          const int blockSize,       /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                          const int numberOfRows,    /* number of rows of the matrix  */
                          const int numberOfColumns, /* number of columns of the matrix */
                          const int* rows,           /* the offset for each row in the vector matrix */
                          const int* columns,        /* the columns number for each matrix entry */
                          const double* matrix,      /* the matrix entries */
                          const double* rhs,         /* the right hand side */ 
                          const int offset)          /* offset for indices (0 for C and C++, 1 for Fortran ) */
{
  XDR xdrs;
  FILE* file = fopen(filename, "wb");

  int bS = blockSize;
  int nR = numberOfRows;
  int nC = numberOfColumns;
  int nZ = rows[ numberOfRows ] - rows[ 0 ];
  int myTest = test;
  char* myMesh = (char *) mesh; 
  char* myName = (char *) nameAndScheme; 

  /* check file */
  if( ! file ) 
  {
    fprintf(stderr,"Couldn't open file `%s' for writing! \n", filename);
    exit( EXIT_FAILURE );
  }

  /* create XDR stream for writing */
  xdrstdio_create( &xdrs, file, XDR_ENCODE);

  /* set xdr pointer */
  setXDRWritePointer( &xdrs );

  /* store dimensions */
  readWriteHeader(readwrite_int_xdr, readwrite_double_xdr, readwrite_string_xdr, 
                  &myName, &myTest, &myMesh,
                  &bS, &nR, &nC, &nZ );

  /* write matrix data to file 
     XDR can only handle non-const data */
  readWriteMatrixAndRhs( readwrite_int_xdr,
                         readwrite_double_xdr,
                         nR, nC,
                         (int *) rows, 
                         (int *) columns, 
                         (double *) matrix, 
                         (double *) rhs, 
                         offset );

  xdrPtr = NULL;

  /* destroy XDR stream */
  xdr_destroy( &xdrs );

  /* close file */
  fclose( file );
}

char* generateFileName(const char* realpath,          /* path to write file to */
                       const char* nameAndScheme, /* filename of file to write to */
                       const int test,            /* number of test , possible set { 1,..., 5 } */
                       const char* mesh)          /* mesh file name */
{
  const char* path = ( realpath == NULL ) ? "." : realpath ;
  const int length = strlen( path ) + strlen( nameAndScheme ) +  strlen( mesh ) + 512 ;
  char* filename = (char *) malloc( length * sizeof( char ) );
  assert( filename );
  
  if( test < 1 || test > 5 ) 
  {
    fprintf(stderr,"ERROR: wrong test number %d , valid are { 1,..., 5} \n", test);  
    exit( 1 );
  }

  sprintf(filename,"%s/%s_test%d_%s.xdr",path,nameAndScheme,test,mesh);
  return filename;
}
                      

/* C/C++ routine */
void writeBenchmarkFile(const char* path,          /* path to write file to */
                        const char* nameAndScheme, /* filename of file to write to */
                        const int test,            /* number of test , possible set { 1,..., 5 } */
                        const char* mesh,          /* mesh file name */
                        const int format,          /* format of CSR storage */
                        const int blockSize,       /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                        const int numberOfRows,    /* number of rows of the matrix  */
                        const int numberOfColumns, /* number of columns of the matrix */
                        const int* rows,           /* the offset for each row in the vector matrix */
                        const int* columns,        /* the columns number for each matrix entry */
                        const double* matrix,      /* the matrix entries */
                        const double* rhs)         /* the right hand side */ 
{
  const int columnoffset = getColumnOffSet( format, numberOfRows, rows, columns );

  char* filename = generateFileName(path,  nameAndScheme, test, mesh );

  doWriteBenchmarkFile(filename, nameAndScheme, test, mesh, blockSize, numberOfRows, numberOfColumns, 
                       rows, columns, matrix, rhs, columnoffset);   

  printf("Wrote file `%s' \n", filename);
  /* free mem */
  free( filename );
}

/* Fortran routine */
void writebenchmarkfile_(const char* path,           /* path to write file to */
                         const char* nameAndScheme,  /* filename of file to write to */
                         const int* test,            /* number of test , possible set { 1,..., 5 } */
                         const char* mesh,           /* mesh file name */
                         const int* format,          /* format of CSR storage */
                         const int* blockSize,       /* block size of DoFs ( = 1 for Finite Volume and Finite Element schemes) */
                         const int* numberOfRows,    /* number of rows of the matrix  */
                         const int* numberOfColumns, /* number of columns of the matrix */
                         const int* rows,            /* the offset for each row in the vector matrix */
                         const int* columns,         /* the columns number for each matrix entry */
                         const double* matrix,       /* the matrix entries */
                         const double* rhs)          /* the right hand side */ 
{
  const int columnoffset = getColumnOffSet( *format, *numberOfRows, rows, columns );

  char* filename = generateFileName(path,  nameAndScheme, *test, mesh );

#if 0
  /* search for spaces */
  const char *search = "  ";
  int fnlength = strlen( fname );
  /* sreach spaces to remove them */
  char* found  = strstr( fname, search );
  char* filename = NULL;

  if ( found )
  {
    fnlength = ((int)( found - fname ));
  }

  filename = (char *) malloc( (fnlength+1) * sizeof(char));
  assert( filename );
  /* copy important parts */
  memcpy( filename, fname, fnlength * sizeof(char));
  /* add terminating character */
  filename[ fnlength ] = '\0';
#endif

  printf("Offsets used: rows %d,  columns %d \n", rows[0], columnoffset );
  /* residual( filename, *numberOfRows, rows, columns, matrix, rhs, solution, offset ); */
  doWriteBenchmarkFile(filename, nameAndScheme, *test, mesh, 
                       *blockSize, *numberOfRows, *numberOfColumns,
                       rows, columns, matrix, rhs, columnoffset);   

  printf("Wrote file `%s' \n", filename);
  /* free mem */
  free( filename );
}

/* read-write solution  */
int readWriteSolution(rw_int_t* rw_int,
                      rw_double_t* rw_double,
                      const int numberOfUnknowns,  /* number of rows of the matrix  */
                      double* solution)            /* the solution */ 
{
  int i, nukwn = numberOfUnknowns;

  rw_int( &nukwn );
  if( nukwn != numberOfUnknowns ) 
  {
    fprintf(stderr,"ERROR: number of unknowns during read of solution wrong! \n"); 
    exit( EXIT_FAILURE );
  }

  /* read-write all offsets of rows */
  for(i=0; i<nukwn; ++i) 
  {
    rw_double( &solution[ i ] );
#ifdef PRINT_OUTPUT
    printf("%d \n", solution[ i ]);
#endif
  }
  return nukwn;
}


void doReadSolution(const char* filename, 
                    const int nUnknowns,
                    double* solution) 
{
  XDR xdrs;
  FILE* file = fopen(filename, "rb");
  int numberOfUnknowns = nUnknowns;

  /* check file */
  if( ! file ) 
  {
    fprintf(stderr,"Couldn't open file `%s' for writing! \n", filename);
    exit( EXIT_FAILURE );
  }

  /* create XDR stream for writing */
  xdrstdio_create( &xdrs, file, XDR_DECODE);

  /* set xdr pointer (read mode) */
  setXDRReadPointer( &xdrs );

  numberOfUnknowns = readWriteSolution(readwrite_int_xdr,
                                       readwrite_double_xdr,
                                       numberOfUnknowns,
                                       solution);

  /* destroy XDR stream */
  xdr_destroy( &xdrs );

  /* close file */
  fclose( file );
}

/* read solution, number of unknowns will be returned */
void readSolution(const char* filename,       /* the file name */
                  const int numberOfUnknowns, /* number of unknowns to be read (for checking) */
                  double* solution)           /* the solution */ 
{
  /* read solution and store number of unknowns */
  doReadSolution( filename, numberOfUnknowns, solution );
}

/* read solution, number of unknowns will be returned */
void readsolution_(const char* fname,           /* the file name */
                   const int* numberOfUnknowns, /* number of unknowns that will be read (stored in the file) */
                   double* solution)            /* the solution */ 
{
  /* search for spaces */
  const char *search = "  ";
  int fnlength = strlen( fname );
  char* found  = strstr( fname, search );
  char* filename = NULL;

  if ( found )
  {
    fnlength = ((int)( found - fname ));
  }

  filename = (char *) malloc( (fnlength+16) * sizeof(char));
  assert( filename );
  /* copy important parts */
  memcpy( filename, fname, fnlength * sizeof(char));
  /* add terminating character */
  filename[ fnlength ] = '\0';

  /* read solution and store number of unknowns */
  doReadSolution( filename, *numberOfUnknowns, solution );

  /* free memory */
  free( filename );
}

/* check residual to make sure structures are ok */
void residual(const char* filename,       /* filename of file */
              const int numberOfRows,     /* number of rows of the matrix  */
              const int* rows,            /* the offset for each row in the vector matrix */
              const int* columns,         /* the columns number for each matrix entry */
              const double* matrix,       /* the matrix entries */
              const double* rhs,          /* the right hand side */ 
              const double* solution,     /* the solution vector */
              const int offset )          /* offset */
{
  /* first entry in rows is also the offset */
  const int rowOffset = rows[ 0 ];
  int i, col;
  double resid = 0.0;
  double res = 0.0;
  for( i = 0; i<numberOfRows; ++i)
  {
    res = -rhs[ i ];
    for( col = rows[ i ] - rowOffset; col < rows[ i+1 ] - rowOffset; ++col )  
    {
      res += matrix[ col ] * solution[ columns[ col  ] - offset ];
    }
    /* printf("res[%d] = %f \n",i,res); */
    resid += (res * res);
  }

  printf("Write benchmark file `%s' with residual %.16e \n", filename, resid);
}

/* read only block size from given file */
int readBlockSize(const char* filename)   /* filename to read blocksize from */
{
  XDR xdrs;
  FILE* file = fopen(filename, "rb");
  int blocksize;

  /* check file */
  if( ! file ) 
  {
    fprintf(stderr,"Couldn't open file `%s' for writing! \n", filename);
    exit( EXIT_FAILURE );
  }

  /* create XDR stream for writing */
  xdrstdio_create( &xdrs, file, XDR_DECODE);

  setXDRReadPointer( &xdrs );

  /* check version of file */
  readWriteVersion( readwrite_double_xdr );

  xdr_int( &xdrs, &blocksize );

  /* destroy XDR stream */
  xdr_destroy( &xdrs );

  /* close file */
  fclose( file );

  return blocksize;
}

/* set xdr pointer to global pointer (write mode) */
void setXDRWritePointer( XDR* xdrs ) 
{
  assert( xdrs != NULL );
  xdrPtr = xdrs ;   
  writeMode = 1 ;
}
void setXDRReadPointer( XDR* xdrs ) 
{
  assert( xdrs != NULL );
  xdrPtr = xdrs ;   
  writeMode = 0 ;
}

#if 0
//**************************************************************************
//**************************************************************************

#include <fstream>
#include <iostream>
#include <stdio.h>

#include <stdlib.h>
#include <math.h>
#include "writexdr.h"

/********************************************************************
 *
 * C program writing example file 
 *
 ********************************************************************/
int read_benchmark_SOL() {
  /* setup system 
     We read the solution of A * x = b with 

          | 2  3  0  0  0 |      | 1 |       |  8 |
          | 3  0  4  0  6 |      | 2 |       | 45 |
      A = | 0 -1 -3  2  0 |  x = | 3 |   b = | -3 |
          | 0  0  1  0  0 |      | 4 |       |  3 |
          | 0  4  2  0  1 |      | 5 |       | 19 |
  */

  /* generated with the call `benchruntime c.xdr' */
  const char* filename = "c.xdr.output.xdr"; 

  const int n = 5 ; /* we know this a-priori since in our programs we 
                       already have the discretization grid */

  double x [5] = { 0.0, 0.0, 0.0, 0.0, 0.0 } ;
  const double u[ 5 ]= { 1.0, 2.0, 3.0, 4.0, 5.0 };

  int i;

  /* call write function, see writexdr.h and writexdr.c for docu 
   * x needs to be pre-allocated with the size to be read */
  readSolution( filename, n, x );

  for( i=0; i<n; ++i )
  {
    if( fabs( x[ i ] - u[ i ] ) > 1e-12 ) 
    {
      fprintf(stderr,"ERROR: solution read (x[ %d ] = %f) is wrong, should be %f !",i,x[i],u[i]);
    }
    printf("x[ %d ] = %f \n",i,x[i]);
  }

  return 0;
}
#endif

/********************************************************************
 *
 * C program writing example file 
 *
 ********************************************************************/
void write_benchmark_CSR( Matrix & MAT, Vector & RHS, string scheme_name, string mesh_name, int test_case ) {

  // setup system
  int offs = 0 ;
  int nrow = MAT.n_rows() ;
  int nnz  = MAT.nnz() ;

  Vector Ax(nnz) ;
  VecInt Ap(nrow+1), Ai(nnz) ;
  
  build_NONSYMM_CSR_format( MAT, Ax, Ap, Ai, offs ) ;
  
  // --------------------------------------------------------------------
  // --------------------------------------------------------------------
  
  const int CSRformat = offs; // default C format = 0
  const int blockSize = 1;    // 1 for most schemes

  const char* path = "."; // path to write file to
  const int   test = test_case ; // number of test, { 1,..., 5 }
  
  const char* name = scheme_name.c_str() ; // name of contributor and scheme
  const char* mesh = mesh_name  .c_str() ; // pass the mesh filename used for calculation
  
  // call write function, see writexdr.h and writexdr.c for docu
  writeBenchmarkFile( path, name, test, mesh,
                      CSRformat, blockSize, nrow, nrow, 
		      Ap.add(), Ai.add(), Ax.add(), RHS.add() ); 
}
