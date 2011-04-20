#define __MINC 40

	module cincom
	    real*8 rpar(__MINC)
	    integer inames
	    integer ipars
	    integer ioldc
	    integer inpar(__MINC)
	    integer iparn(__MINC)
	    integer inapa(__MINC)
	    integer iargs
	    integer ipmls
	    integer iolbuf
	    integer lstpar
	    integer lstnam
	end module cincom

	module cincoc
		character*8 comand
		character*8 vname(__MINC)
		character*1024 title
		character*1024 reslin
		character*1024 inline
		character*20 arglst(__MINC)
		character*20 pmlist(__MINC,2)
		character*1024 rlbuf
	end module cincoc

	module icpathes
		character*1024 data_path
		character*1024 save_path
		character*1024 makro_path
		character*1024 home
		character*1024 datreat_path,PWD
		character*1024 history(0:20)
	end module icpathes

	module constants
		integer, parameter :: minc=__MINC
	end module constants
