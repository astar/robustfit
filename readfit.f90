MODULE fits
CONTAINS
SUBROUTINE readfit(filename, wave, flux)
integer :: status, blocksize, naxis, pcount,gcount, bitpix, nrows
integer :: group = 1, extver = 0, frow = 1, felem = 1, hdutype
real :: nullval = 0.0
logical :: simple,extend,anyf
real, dimension(:), allocatable :: wave,flux
CHARACTER(len=*) :: filename
status = 0

!call get_command_argument(1, filename)

call ftopen(15,filename,1,blocksize,status)
call FTMAHD(15,2, hdutype,status)
call ftgnrw(15,nrows,status)

allocate(wave(nrows),flux(nrows))


call ftgcve(15,1,frow,felem,size(wave),nullval,wave,anyf,status)
call ftgcve(15,2,frow,felem,size(flux),nullval,flux,anyf,status)
END SUBROUTINE readfit
END MODULE fits
