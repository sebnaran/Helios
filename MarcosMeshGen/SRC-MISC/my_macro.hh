# ifndef _MY_MACRO_HH
# define _MY_MACRO_HH

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// INPUT/OUTPUT Macros
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

// to get namespace std
# include<iostream> 
# define OUTPUT std::cout
# define INPUT  std::cin
# define ERRPUT std::cerr

# define PRT(A)      ( OUTPUT << (#A) << " = " << (A) << std::endl )
# define VAL(A)      ( OUTPUT << (#A) << " = " << (A) << " " << std::flush )
# define VALV(I,A)   ( OUTPUT << "  " << (#A) << "(" << I << ") =" << A((I)) << " " <<std::flush )
# define VALA(I,A)   ( OUTPUT << "  " << (#A) << "[" << I << "] =" << A[(I)] << " " <<std::flush )
# define VALM(I,J,A) ( OUTPUT << "  " << (#A) << "(" << I << "," << J << ") =" << A((I),(J)) << std::flush )
# define PRTA(I,A)   ( OUTPUT << "  " << (#A) << "[" << I << "] =" << ( abs(A[(I)])<1.e-20? 0. : A[(I)] )<< std::endl )
# define PRTV(I,A)   ( OUTPUT << "  " << (#A) << "(" << I << ") =" << ( abs(A((I)))<1.e-20? 0. : A((I)) )<< std::endl )
# define PRTM(I,J,A) ( OUTPUT << "  " << (#A) << "(" << I << "," << J << ") =" << ( abs(A((I),(J)))<1.e-20? 0. : A((I),(J)) )<< std::endl )

# define DBGF(A)     do { OUTPUT << #A << std::endl << std::flush ; } while(0)
# define MSG(A)      do { OUTPUT << A ; } while(0)
# define MSGF(A)     do { OUTPUT << A << std::flush << std::endl ; } while(0)
# define INSERT(A)   do { OUTPUT << "insert " << (#A) << " = " ; INPUT >> A ; } while(0)

const bool _verbose_ = bool(false) ;
# define VRBS(A)     do { if ( _verbose_ ) { OUTPUT << #A << std::flush << std::endl ; } } while(0)

# define PAUSE								\
  do {									\
    char dummy=' ' ;							\
    ERRPUT << "pause (0 to continue)" ;					\
    while ( dummy != '0' ) { INPUT >> dummy ; }				\
  } while(0)

# define indexType int
# define PRT_VEC(V)     for( indexType i=0; i<V.size(); ++i ) { PRTV(i,V) ; }
# define PRT_ARR(V)     for( indexType i=0; i<V.size(); ++i ) { PRTA(i,V) ; }
# define PRT_ARRAY(N,V) for( indexType i=0; i<N; ++i )        { PRTA(i,V) ; }
# define PRT_MATRIX(MAT)				     \
  do {							     \
    for ( indexType j=0 ; j<MAT.size_j() ; ++j ) {	     \
      MSG("print col --> ") ; PRT(j) ;			     \
      for ( indexType i=0 ; i<MAT.size_i() ; ++i ) {	     \
	PRTM(i,j,MAT) ;					     \
      }							     \
    }							     \
  } while(0)

# define PRT_MATROW(MAT)				     \
  do {							     \
    for ( indexType i=0 ; i<MAT.size_i() ; ++i ) {	     \
      MSG("print row --> ") ; PRT(i) ;			     \
      for ( indexType j=0 ; j<MAT.size_j() ; ++j ) {	     \
	PRTM(i,j,MAT) ;					     \
      }							     \
    }							     \
  } while(0)

# define LINE(A)							\
  do {									\
    for( indexType i=0; i<20; ++i ) { OUTPUT << (#A) ; }		\
    OUTPUT << std::endl ;						\
  } while(0)

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template<class T>
void swap_items( T & a, T & b ) {
  T tmp = a ;
  a = b ;
  b = tmp ;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#if 0
namespace eat {

  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }

  static istream & eatchar(istream & s) { s.get() ; return s ; }

  static istream & eatcomments(istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
            && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }

  void fatal_error(string str) {
    cerr << "fatal error\n" ;
    cerr << str << endl << flush ;
    exit(0) ;
  }

}
using namespace eat ;
#endif

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Timing stuff
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#ifndef TIMER_H_DEF
#define TIMER_H_DEF

// Unix based system specific
#include <sys/time.h>

class Timer {
public:
  Timer() ;                                    // default constructor
  ~Timer() ;                                   // default destructor
  
  void   start() ;                             // start timer
  void   stop () ;                             // stop the timer
  double getElapsedTime() ;                    // get elapsed time in second
  double getElapsedTimeInSec() ;               // get elapsed time in second (same as getElapsedTime)
  double getElapsedTimeInMilliSec() ;          // get elapsed time in milli-second
  double getElapsedTimeInMicroSec() ;          // get elapsed time in micro-second
  
  void report( const char * what ) ;

private:
  double  startTimeInMicroSec ;                // starting time in micro-second
  double  endTimeInMicroSec ;                  // ending time in micro-second
  int     stopped ;                            // stop flag 
  timeval startCount ;                         // ticks per second
  timeval endCount ;                           //
};

Timer::Timer() {
  startCount.tv_sec = startCount.tv_usec = 0 ;
  endCount.tv_sec   = endCount.tv_usec = 0 ;
  stopped = 0 ;
  startTimeInMicroSec = 0 ;
  endTimeInMicroSec   = 0 ;
}

Timer::~Timer() {}

// start timer.
// startCount will be set at this point.
void Timer::start() {
  stopped = 0 ; // reset stop flag
  gettimeofday( &startCount, NULL ) ;
}

// stop the timer.
// endCount will be set at this point.
void Timer::stop() {
  stopped = 1; // set timer stopped flag
  gettimeofday( &endCount, NULL ) ;
}

// compute elapsed time in micro-second resolution.
// other getElapsedTime will call this first, then convert to correspond resolution.
double Timer::getElapsedTimeInMicroSec() {
  if(!stopped) {
    gettimeofday(&endCount, NULL) ;
  }
  
  startTimeInMicroSec = ( startCount.tv_sec * 1000000.0 ) + startCount.tv_usec ;
  endTimeInMicroSec   = ( endCount.tv_sec   * 1000000.0 ) + endCount.tv_usec ;

  return endTimeInMicroSec - startTimeInMicroSec;
}

// divide elapsedTimeInMicroSec by 1000
double Timer::getElapsedTimeInMilliSec() {
  return this->getElapsedTimeInMicroSec() * 0.001 ;
}

// divide elapsedTimeInMicroSec by 1000000
double Timer::getElapsedTimeInSec() {
  return this->getElapsedTimeInMicroSec() * 0.000001 ;
}

// same as getElapsedTimeInSec()
double Timer::getElapsedTime() {
  return this->getElapsedTimeInSec() ;
}

void Timer::report( const char * what ) { 
  OUTPUT << "Timing for " << what <<" = "
	 << this->getElapsedTimeInMicroSec() << " micro seconds\n" ;
}

# define TIMING(A) {{							\
      OUTPUT << "running " << #A << std::endl << std::flush ;		\
      Timer o_clock ;							\
      o_clock . start() ;						\
      (A) ;								\
      o_clock . stop() ;						\
      OUTPUT << "end of  " << #A  << std::endl << std::flush ;		\
      o_clock . report( #A ) ; 						\
      LINE(---) ;                                                       \
    }}

#endif // TIMER_H_DEF

// print time prints the current date in YMDHMS format
void print_time() {
  const int TIME_SIZE=40 ;

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// to convert integers to strings
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#if 0
// C-like implementation !!!
# ifndef __FOO__
# define __FOO__
# include<stdio.h>
void foo (char buf[], unsigned n) {
  sprintf(buf,"%d",n);
  puts(buf); // will print n...
}

# endif // end of __FOO__
#endif

#include<sstream>
string append_int_to_string( int n, string inp_str ) {
  stringstream outp_stream(inp_str.c_str()) ;
  outp_stream.seekp(0, std::ios::end);
  outp_stream << n ;
  return outp_stream.str() ;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

// usage:
// ======
// vector<sort_pair> list_pair ;
// list_pair.push_back( sort_pair(i0,v0) ) ;
// sort( list_pair.begin(), list_pair.end() ) ;

struct sort_pair {
public:
  static double eps ;
  int i ;
  double v ;
  sort_pair( int _i, double _v ) : i(_i), v(_v) {}
  ~sort_pair() {}
  static double _abs( double x ) { return x>0. ? x : -x ; }
} ;

double sort_pair::eps = 1.e-12 ;

bool operator< ( const sort_pair & p0, const sort_pair & p1 ) { return p0.v<p1.v ; }
bool operator==( const sort_pair & p0, const sort_pair & p1 ) { return p0.i==p1.i && sort_pair::_abs(p0.v-p1.v)<sort_pair::eps ; }
bool operator!=( const sort_pair & p0, const sort_pair & p1 ) { return !( p0==p1 ) ; }


# endif // end of _MY_MACRO_HH
