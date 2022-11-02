program fourier

  real, parameter :: pi = 3.141592654
  
  logical :: ex
  
  character(len=30) :: get(20), signal

  integer :: n_index, i, j, ierr

  real :: aux, aux2, n
  real, dimension(1000000) ::func, re, im, ref, imf


  write(*, *) "  .--.--.                   ,---,                      ,---,. "
  write(*, *) " /  /    '.                '  .' \            ,---.  ,'  .' | "
  write(*, *) "|  :  /`. /          ,--, /  ;    '.         /__./|,---.'   | "
  write(*, *) ";  |  |--`         ,'_ /|:  :       \   ,---.;  ; ||   |   .' "
  write(*, *) "|  :  ;_      .--. |  | ::  |   /\   \ /___/ \  | |:   :  |-, "
  write(*, *) " \  \    `. ,'_ /| :  . ||  :  ' ;.   :\   ;  \ ' |:   |  ;/| "
  write(*, *) "  `----.   \|  ' | |  . .|  |  ;/  \   \\   \  \: ||   :   .' "
  write(*, *) "  __ \  \  ||  | ' |  | |'  :  | \  \ ,' ;   \  ' .|   |  |-, "
  write(*, *) " /  /`--'  /:  | : ;  ; ||  |  '  '--'    \   \   ''   :  ;/| "
  write(*, *) "'--'.     / '  :  `--'   \  :  :           \   `  ;|   |    \ "
  write(*, *) "  `--'---'  :  ,      .-./  | ,'            :   \ ||   :   .' "
  write(*, *) "             `--`----'   `--''               '---: |   | ,'   "
  write(*, *) "                                                   `----'     "
  write(*, *) ""
  write(*, *) ""  
  write(*, *) "                       ** s_filter **"
  write(*, *) ""
  write(*, *) "             Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent" 
  write(*, *) "Properties in Chemical Interfaces. Denys E. S. Santos," 
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho," 
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model. 2019."
  write(*, *) ""
  write(*, *) ""
  write(*, *) "s_filter performes a Discrete Fourier Transform (DFT)"
  write(*, *) "on the signal described in the input file, so that the"
  write(*, *) "user becomes able to filter noise inside the data"
  write(*, *) ""
  write(*, *) "Usage: s_filter -in file.dat"
  write(*, *) ""
  write(*, *) "file.dat ---- contains the data to be filtered"

  signal = 'miss'

  do i=1, 20
     
     call getarg(i, get(i))
     
  end do

  do i=1, 20
     
     if (get(i)=='-in')then
        
        signal = get(i+1)

     end if
     
  end do

  if (signal=='miss')then

     write(*, *)
     write(*, *)'Input file is missing'
     write(*, *)
     stop

  end if

  open(1, file=signal, status='old', iostat=ierr)

  if(ierr /= 0) then
     
     write(*, *)
     write(*, *) ' Unable to open file ', signal
     write(*, *)
     stop
     
  endif

  inquire(file='transf.xvg', exist=ex)

  if (ex) then
     
     call execute_command_line("mv transf.xvg transf_p.xvg")
     write(*, *)
     write(*, *) ' Previous transf.xvg file is now backed up to transf_p.xvg !'
     
  end if
  
  open(2, file='transf.xvg', status='new', iostat=ierr)
  
  if(ierr /= 0) then
     
     write(*, *)
     write(*, *) ' Unable to open file transf.xvg'
     write(*, *)
     stop
     
  endif
  
  write(2, *) '@    title "Fourier Transform"'
  write(2, *) '@    xaxis  label "Hz"'
  write(2, *) '@    yaxis  label "|F[n]|"'

  
  inquire(file='filtered.xvg', exist=ex)

  if (ex) then

     call execute_command_line("mv filtered.xvg filtered_p.xvg")
     write(*, *)
     write(*, *) ' Previous filtered.xvg file is now backed up to filtered_p.xvg !'

  end if

  open(3, file='filtered.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) ' Unable to open file filtered.xvg'
     write(*, *)
     stop

  endif

  write(3, *) '@    title "Inverse Fourier Transform"'
  write(3, *) '@    xaxis  label "t"'
  write(3, *) '@    yaxis  label "|f[t]|"'
  
  do i=1, 1000000

     func(i) = 0
     re(i) = 0
     im(i) = 0

  end do
  
  !===================================

  n_index = 1
  
  do while (ierr>=0)
     
     read(1, *, iostat=ierr) aux, aux2
     
     if (ierr == 0)then

        func(n_index) = aux2
        n_index = n_index + 1

     end if
     
  end do

  n_index = n_index - 1

  write(*, *)
  write(*, *) 'Maximum frequency = ', n_index, '1/cicle'
  write(*, *)

  if(n_index>1000000) then

     write(*, *)
     write(*, *) ' Too many points'
     write(*, *)
     stop

  end if
  
  
  !===================================
  ! calculando a transformada
  
  do i=0, n_index -1

     re(i+1) = 0
     im(i+1) = 0
     
     do j=0, n_index -1 

        re(i+1) = re(i+1) + func(j+1)*cos(2*pi*i*j/n_index)
        im(i+1) = im(i+1) - func(j+1)*sin(2*pi*i*j/n_index)

     end do

     aux = i

        write(2, *, iostat=ierr) aux, sqrt(re(i+1)*re(i+1) + im(i+1)*im(i+1))

     if(ierr>0) then

        write(*, *)
        write(*, *) ' Problem by writing output '
        write(*, *)
        stop

     end if
     
  end do

  close(2)

  !==inserindo brancos para frequencias altas===
  do i=int(n_index/4 +1), n_index

     re(i) = 0
     im(i) = 0

  end do

  !==============================================
  ! calculando a funcao atraves da transformada

  do j=0, n_index - 1

     ref(j+1) = 0
     imf(j+1) = 0

     do i=0, n_index - 1

        ref(j+1) = ref(j+1) + re(i+1)*cos(2*pi*i*j/n_index) - im(i+1)*sin(2*pi*i*j/n_index)
        imf(j+1) = imf(j+1) + re(i+1)*sin(2*pi*i*j/n_index) + im(i+1)*cos(2*pi*i*j/n_index)

     end do

     aux = j

!    write(3, *, iostat=ierr) aux, sign(1.00 ,atan(imf(j+1)/ref(j+1)))*sqrt(ref(j+1)*ref(j+1) + imf(j+1)*imf(j+1))/n_index
     write(3, *, iostat=ierr) aux, sqrt(ref(j+1)*ref(j+1) + imf(j+1)*imf(j+1))/n_index                    

     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by writing output '
        write(*, *)
        stop

     end if

  end do
  
  close(1)
  close(3)

  write(*, *)
  write(*, *) " Finished"
  write(*, *)
    
end program fourier
