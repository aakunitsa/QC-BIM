module consts
  
  implicit none
  
  real*8, parameter :: boltz    = 0.316682968d-5 ! E_h/K
  real*8, parameter :: bohr2ang = 0.529177249d0  ! A/bohr
  real*8, parameter :: ang2bohr = 1.0d0/bohr2ang ! bohr/A: A --> bohr
  real*8, parameter :: au_time  = 0.024188843d0  ! fs/A.U. 
  real*8, parameter :: amu2au   = 1822.8885300625575d0 ! at.weigh --> A.U.
  real*8, parameter :: au_bar   = 2.942191d8  ! a.u.(E_h/bohr**3) --> bar
  !real*8, parameter :: HBAR =0.0151704885d0 ! kcal/mol * ps
  
  !,,,,,,
  real*8, parameter ::     PI= 3.14159265358979323844d0
  real*8, parameter ::  TWOPI= 6.28318530717958647688d0
  real*8, parameter :: FOURPI=12.56637061435917295376d0
  real*8, parameter ::    TSP= 1.12837916709551257390d0 !  2/sqrt(pi)
  real*8, parameter ::   SQPI= 1.77245385090551602729d0 ! sqrt(pi)
  real*8, parameter :: RADIAN= 57.29577951308232087721d0
  
end module consts

