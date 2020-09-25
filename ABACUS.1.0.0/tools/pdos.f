	program fmPDOS
C
C       extracts data from PDOS file of SIESTA,
C       taking care of flexible format in index, atom_index, etc. fields
C       Andrei Postnikov, University Paul Verlaine - Metz
C       Jul 2010
C       postnikov@univ-metz.fr
C
      implicit none
      integer ii1,io1
      parameter (ii1=11,io1=12)
      integer nt,nmax,i0,ispin,nspin,norbs,it,npts,nene,is,idos,l,lref,
     .        m,mode,mref,z,atind,indref,index,nline,
     .        iquoted,iparsed
      parameter (nmax=10000)
      double precision ene(nmax),dos(nmax,4),dos1(4)
      character inpfil*60,outfil*60,string*80,llabel*80,rlabel*80,
     .          species*6,squoted*6,chlab*6,owrite*1,s1*1
      logical filexist,redos,nptdef
      external iquoted,squoted,iparsed

      do it=1,nmax
      do is=1,2
        dos(it,is)=0.d0
      enddo 
      enddo 
      nptdef = .false.     ! We don't know whether the file contains <npoints>

C --- open input and outpute files: ----------------------
      write (6,*) ' Input file name (PDOS):'
      read (5,*) inpfil
      inquire (file=trim(inpfil), exist=filexist)
      if (filexist) then
        open (ii1,file=trim(inpfil),status='old',form='formatted')
      else
        write (6,*) ' File does not exist!'
        stop
      endif

      write (6,*) ' Output file name :'
      read (5,*) outfil
      inquire (file=trim(outfil), exist=filexist)
      if (filexist) then
        write (6,*) ' File ',trim(outfil),' exists. Overwrite? (Y/N)'
        read (5,*) owrite
        if (owrite.eq.'Y'.or.owrite.eq.'y') then
          open (io1,file=trim(outfil),form='formatted',status='REPLACE')
        else
          write (6,*) ' Then rename is first. Bye...'
          stop
        endif
      else
        open (io1,file=trim(outfil),form='formatted',status='NEW')
      endif

C --- select projected DOS which are needed: -------------
  100 continue
      write (6,*) ' Extract data for atom index',
     .            ' (enter atom NUMBER, or 0 to select all),' 
      write (6,*) ' or for all atoms of given species',
     .            ' (enter its chemical LABEL):'
      read (5,*) string
C     check first character of string, whether it is decimal
      i0 = ichar(string(1:1))
      if (i0.ge.48.and.i0.le.57) then
C       1st character passed is a number;
C       assume the whole string is atom number
        mode = 1
        read (string,*,err=401) indref
        write (io1,301)indref
      elseif ((i0.ge.65.and.i0.le.90).or.(i0.ge.97.and.i0.le.122)) then
C       1st character passed is a character [A-Z,a-z];
C       assume the whole string is atom label
        mode = 2
        read (string,*,err=402) species
        write (io1,302) trim(species)
      else
C       1st character is a special symbol; presumably an error
        write (6,306) i0
        goto 100
      endif

              write (6,*) ' Extract data for l= ... (-1 for all l ):'
        read (5,*) lref
        if (lref.ne.-1) then 
          write (6,*) ' Extract data for m= ... (9 for all m ):'
          read (5,*) mref
        else
          mref=9
        endif 
    
      redos = .false.
      nene=0
      nline = 0

C --- line by line read and analyze PDOS file: -----------
   11 continue
      nline = nline+1
      read (ii1,'(a80)',err=201,end=202) string
      if (string(1:6).eq.'<pdos>') goto 11
      if (string(1:7).eq.'<nspin>') then
        llabel = '<nspin>'
        rlabel = '</nspin>'
        nspin = iparsed(string,llabel,rlabel)
        if (nspin.eq.8) then
          write (6,*) ' nspin=8 changed to nspin=4'
          nspin=4
        endif
        goto 11
      elseif (string(1:11).eq.'<norbitals>') then
        llabel = '<norbitals>'
        rlabel = '</norbitals>'
        norbs = iparsed(string,llabel,rlabel)
        goto 11
      elseif (string(1:9).eq.'<npoints>') then
        llabel = '<npoints>'
        rlabel = '</npoints>'
        nptdef = .true.
        npts = iparsed(string,llabel,rlabel)
        goto 11
      elseif (string(1:8).eq.'<energy_') then   !  list of energies opens:
        nene=0
        goto 11
      elseif (string(1:8).eq.'        ') then   !  must be an energy line:
                                                !  read energy values
        nene=nene+1
        read (string,*) ene(nene)
        goto 11
      elseif (string(1:16).eq.'</energy_values>') then ! list of energies closes
        nt=nene
C       write (6,*) nt,'  energy values found:'
C       do it=1,nt
C         write(6,203) it,ene(it)
C 203     format(i5, f12.6)
C       enddo
        if (nt.gt.nmax) then
          write(6,*)'  nt=',nt,' > nmax=',nmax
          stop
        endif 
        if (nptdef.and.nt.ne.npts) then
          write(6,*)'  nt=',nt,' differs from npoints=',npts
          stop
        endif 
        goto 11

      elseif (string(1:8).eq.'<orbital') then      !  new orbital follows:
        goto 11
      elseif (string(2:7).eq.'index=') then        !  orbital index line:
        index = iquoted(string,nline)
        goto 11
      elseif (string(2:12).eq.'atom_index=') then  !  atom index line:
        atind = iquoted(string,nline)
        if (mode.eq.1.and.(atind.eq.indref.or.indref.eq.0)) redos=.true.
        goto 11
      elseif (string(2:9).eq.'species=') then      !  atom species line:
        chlab = squoted(string,nline)
C       write (6,*) ' Compare chlab=',trim(chlab),
C    .              '-- and specias=',trim(species),'---'
        if (mode.eq.2.and.(trim(chlab).eq.trim(species))) then
          redos=.true.
C         write (6,*) ' They match'
        endif
        goto 11
                 elseif (string(2:3).eq.'l=') then   !  l-value:
        l = iquoted(string,nline)
        if (lref.ne.-1.and.l.ne.lref) redos=.false.
        goto 11
      elseif (string(2:3).eq.'m=') then   !  m-value:
        m = iquoted(string,nline)
        if (mref.ne.9.and.m.ne.mref) redos=.false.
        goto 11
      elseif (string(2:3).eq.'z=') then   !  z-value:
        z = iquoted(string,nline)
        read (ii1,'(a1)') s1    !   read in an extra line (closing > )
        goto 11
      elseif (string(1:6).eq.'<data>') then        !  list of PDOS follows:
        if (redos) write (io1,303) atind,l,m,z
        do it=1,nt
          read (ii1,*) (dos1(ispin),ispin=1,nspin)
          if (redos) then
            do ispin=1,nspin
              dos(it,ispin)=dos(it,ispin)+dos1(ispin)
            enddo
          endif
        enddo
        redos=.false.
        goto 11
      elseif (string(1:7).eq.'</data>') then       !  list of PDOS closes:
        goto 11
      elseif (string(1:10).eq.'</orbital>') then   !  orbital ends:
        goto 11
      elseif (string(1:7).eq.'</pdos>') then       !  regular end
        goto 203
      else
        write(6,*) ' Unknown identifier in PDOS file, line ',nline
        write(6,*) string
        stop
      endif
     
  203 continue
      close (ii1)

C --- write down accumulated DOS values: -----------------
      write (io1,304)
      do it=1,nt
        write (io1,305) ene(it),(dos(it,ispin),ispin=1,nspin)
      enddo
      close (io1)
      stop

  201 continue
      write(6,*) 'Error reading PDOS file, line ',nline
      stop
  202 continue
      write(6,*) 'Unexpected end of PDOS file, line ',nline
      stop
  401 write (6,*) ' Error reading ',trim(string),' as numeric'
      goto 100
  402 write (6,*) ' Error reading ',trim(string),' as string'
      goto 100
  301 format('#',/,'#   partial DOS for atom index ',i5,/,'#')
  302 format('#',/,'#   partial DOS for atom species: ',(a),/,'#')
  303 format('#   Add data for atom_index =',i4,',  n,l,m,z=',4i3)
  304 format('#',/,'#    Energy',10x,'spin 1',8x,'spin2',/,'#')    
  305 format(f13.7,4f15.8)
  306 format(' Illegal first character (ASCII =',i4,').',/
     .       ' Atom number must start from [0-9],', 
     .       ' atom label - from [A-Z,a-z]. Try again...')
      end
C
C ......................................................................
C
      integer function iquoted(string,nline)
C
C     returns integer value contained in string between two delimiters
C
      implicit none
      integer lquote,rquote,nline 
      character*80 string
      character*1  delim
      data delim /'"'/
      lquote = scan(string,delim)  
      rquote = scan(string,delim,back=.true.)
      if (lquote.lt.rquote) then
        read (string(lquote+1:rquote-1),*) iquoted
      else               
        write (6,*) ' Error locating quotes in line ',nline
        write (6,*) ' lquote, rquote =',lquote, rquote
        write (6,*) string
        stop
      endif              
      return
      end
C
C ......................................................................
C
      character*6 function squoted(string,nline)
C
C     returns substring contained in string between two delimiters
C
      implicit none
      integer lquote,rquote,nline 
      character*80 string
      character*1  delim
      data delim /'"'/
      lquote = scan(string,delim)  
      rquote = scan(string,delim,back=.true.)
      if (lquote.lt.rquote) then
        read (string(lquote+1:rquote-1),*) squoted
      else               
        write (6,*) ' Error locating quotes in line ',nline
        write (6,*) ' lquote, rquote =',lquote, rquote
        write (6,*) string
        stop
      endif              
      return
      end
C
C ......................................................................
C
      integer function iparsed(string,llabel,rlabel)
C
C     reads in an integer value from the field of string 'string'
C     which is situated between two substrings 'llabel' and 'rlabel'.
      implicit none
      integer first,last 
      character*80 string
      character*80 llabel,rlabel
      
C --- position in 'string' just after appearance of substring 'llabel'
      first = scan(string,llabel(1:len_trim(llabel))) +
     +        len_trim(llabel) 
C --- position in 'string' just before the appearance of substring 'rlabel'
      last = scan(string,rlabel(1:len_trim(rlabel)),back=.true.) -
     -       len_trim(rlabel)
      if (first.le.last) then
        read (string(first:last),*,err=801) iparsed
        return
      else               
        write (6,*) ' Fail to locate integer field in the string :'
        write (6,*)  string
        write (6,*) ' first, last =',first,last
        stop
      endif              
  801 write (6,*) ' Fail to read integer number from positions ',
     .             first,'  through ',last,' of the string :'
      write (6,*) string
      stop
      end
	
