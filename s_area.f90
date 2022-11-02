module area

  implicit none

  real, parameter :: pi = 3.141592654
  
  type vet1

     character(len=50) :: atom, resid, ident, code
     real ::  x, y, z
     integer :: n_atom, n_resid

  end type vet1

  type vet2

     real :: x, y, z

  end type vet2


end module area

program ar_grid
  
  use area

  logical :: ex, bin, outer, p_grid, rmsd, l_coarse, begin, end, skip
  logical :: eval_skip

  integer :: i, j, ierr, n_grid, frame, aux, bin_out, noi1, noi2
  integer :: n_index, n_lipid, num, num2, i_atom, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip
  integer, dimension(50000) :: in_num, in_num2

  real :: dist_x, dist_y, dist_z, s_soma, hour, minu, sec, noir
  real :: x_max, x_min, y_max, y_min, dx, dy, dist, s_grid, al
  real :: s_area, la, lb, lc, start, finish, r_fit, aux2, gridx, gridy

  character(len=30) :: coord,  ind, atom, ind2
  character(len=30), dimension(30) :: get

  type(vet1), dimension(50000) :: store, store2, coarse, coarse2
  type(vet1) :: buff
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat1, mat2
  type(vet2), dimension(1001,1001) :: grid, grid2, grid3

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
  write(*, *) "                       ** s_area **"
  write(*, *) ""
  write(*, *) "             Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent" 
  write(*, *) "Properties in Chemical Interfaces. Denys E. S. Santos," 
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho," 
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model. 2019."
  write(*, *) ""
  write(*, *) ""
  write(*, *) "s_area calculates the time-dependent area per molecule of the" 
  write(*, *) "interface."
  write(*, *) ""
  write(*, *) "Usage: s_area -in file.pdb -ind1 file1.ndx -ind2 file2.ndx"
  write(*, *) ""
  write(*, *) "file.pdb ---- atomic coordinates in PDB format"
  write(*, *) "file.ndx ---- index file containing user-selected atoms used" 
  write(*, *) "to fit the grid points to the chemical surface."
  write(*, *) ""
  write(*, *) "Options:"
  write(*, *) ""
  write(*, *) "-bin             defines the number of rectangular partition "
  write(*, *) "bins along the x- and y-axes"
  write(*, *) ""
  write(*, *) "-outer           automatically selects the "
  write(*, *) "surface/interface outmost atoms to fit the grid points. This" 
  write(*, *) "option overwrites the user-selected index files."
  write(*, *) ""
  write(*, *) "-grid            generates a PDB file containing the grid" 
  write(*, *) "points used in the fitting for the last frame in the "
  write(*, *) "trajectory file."
  write(*, *) ""
  write(*, *) "-rmsd            calculates the RMSD between the fitted" 
  write(*, *) "grid and the selected atoms in the index files. This" 
  write(*, *) "estimates how precisely is the grid surface fitted to the" 
  write(*, *) "chemical surface throughout the trajectory file."
  write(*, *) ""
  write(*, *) "-coarse          generates a coarse grid over the surface" 
  write(*, *) "index atoms from which a finer grid will be generated. This" 
  write(*, *) "recommended for surfaces defined by atoms which greatly" 
  write(*, *) "fluctuate throughout the trajectory. "
  write(*, *) ""
  write(*, *) "-begin           first frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-end             last frame to use in the calculations"
  write(*, *) ""
  write(*, *) "-skip            number of trajectory frames to be skipped" 
  write(*, *) "during the analysis "
  
!
! pegando os arquivos de entrada
!

  outer = .false.
  bin = .false.
  p_grid = .false.
  coord = 'miss'
  ind = 'miss'
  ind2 = 'miss'
  rmsd = .false.
  l_coarse = .false.
  begin = .false.
  end = .false.
  skip = .false.
  
  do i=1, 20

     call getarg(i, get(i))

  end do

  do i=1, 20

     if (get(i)=='-in')then

        coord = get(i+1)

     end if

     if (get(i)=='-ind1')then

        ind = get(i+1)

     end if

     if (get(i)=='-ind2')then

        ind2 = get(i+1)

     end if
     
     if (get(i)=='-bin')then

        bin= .true.

     end if

     if (get(i)=='-outer')then

        outer= .true.

     end if

     if (get(i)=='-grid')then

        p_grid = .true.

     end if

     if (get(i)=='-rmsd')then

        rmsd = .true.

     end if
     
     if (get(i)=='-coarse')then

        l_coarse = .true.

     end if

     if (get(i)=='-begin')then

        begin = .true.
                
     end if

     if (get(i)=='-end')then

        end = .true.

     end if

     if (get(i)=='-skip')then

        skip = .true.

     end if
     
  end do
  
  if (coord=='miss')then

     write(*, *)
     write(*, *)'PDB file is missing'
     write(*, *)
     stop

  end if

  if (.not.outer)then

     if (ind=='miss')then

        write(*, *)
        write(*, *)'First index file is missing'
        write(*, *)
        stop

     end if

     if (ind2=='miss')then

        write(*, *)
        write(*, *)'Second index file is missing'
        write(*, *)
        stop

     end if

  end if
  
  !==Lendo os arquivos de index=========================
  if (.not.outer)then
     
     open(2, file=ind, status='old', iostat=ierr)
     
     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file ', ind
        write(*, *)
        stop
        
     endif
     
     open(3, file=ind2, status='old', iostat=ierr)
     
     if(ierr/=0) then

        write(*, *)
        write(*, *) 'Unable to open file', ind2
        write(*, *)
        stop
        
     end if

     ! Inicio da leitura do index1====================

     ierr = 0
     n_index = 1
     
     do while (ierr==0)
        
        read(2, *, iostat=ierr) aux
        
        if (ierr==0) then
           
           in_num(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', ind
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num(i) = -1
        
     end do

     write(*, *)
     write(*, *) ind, " input file is OK"
     
     ! Fim da leitura do index ==========================

     close(2)
     
     ! Inicio da leitura do ind2 ========================

     ierr = 0
     n_index = 1
     
     do while (ierr==0)
        
        read(3, *, iostat=ierr) aux
        
        if (ierr==0)then
           
           in_num2(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', ind2
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num2(i) = -1
        
     end do

     write(*, *)
     write(*, *) ind2, " input file is OK"
     
     ! Fim da leitura do ind2===================

     close(3)

  end if
  
2 format(a10)
  write(*, *)
  write(*, 2, advance='no') ' lipid  = '
  read(*, *) n_lipid
  
  if ((bin).or.(outer)) then

     write(*, *)
     write(*, 2, advance='no') '    bin = '
     read(*, *) n_grid

  end if

  if (outer) then

     write(*, *)
     write(*, 2, advance='no') 'bin_out = '
     read(*, *, iostat=ierr) bin_out

  end if

  !=================definindo frames para inicio e fim========
  !=================definindo skip============================
  
  frame = 0
  
  if (begin) then

     write(*, *)
     write(*, 2, advance='no') " Begin =  "
     read(*, *) fr_in

  else

     fr_in = 1

  end if
  
  if (end) then

     write(*, *)
     write(*, 2, advance='no') " End   =  "
     read(*, *) fr_end

  else

     fr_end = fr_in + 2 !==para que fr_end seja maior que frame

  end if

  if (skip) then

     write(*, *)
     write(*, 2, advance='no') " n_skip = "
     read(*, *) n_skip

  else

     n_skip = 0
     
  end if
  !============================================================
  
  !============================================================
  ! calculando o espaçamento

  if (.not.outer)then

     bin_coarse = nint(sqrt(real(n_index-1)) - 1)
     
     if (.not.bin)then
        
        n_grid = nint(sqrt(real(n_index-1)) - 1)
        write(*, *)
        write(*, *) 'STD_BIN = ', n_grid
        
     end if
     
  else

     bin_coarse = bin_out
     
  end if
  !===================================================
  
  !===================================================
  ! definição da matriz de pontos externos
  
  if (outer)then
     
     allocate (mat1(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'
     
     allocate (mat2(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'

     mat1(1:bin_out+1,1:bin_out+1)%z = -1000
     mat2(1:bin_out+1,1:bin_out+1)%z = 1000
          
  end if
  
  !=================================================

  write(*, *)
  write(*, *) 'Processing ..........' !==============================iniciando processamento

  open(1, file=coord, status='old', iostat=ierr)
  
  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file',coord
     write(*, *)
     stop
     
  endif
  
  inquire(file='area-med.xvg', exist=ex)
  
  if (ex) then
     
     call execute_command_line("mv area-med.xvg area-med_p.xvg")
     write(*, *)
     write(*, *) 'Previous area-med.xvg file is now backed up to area-med_p.xvg !'

  end if

  open(2, file='area-med.xvg', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file area-med.xvg'
     write(*, *)
     stop

  endif

  write(2, *) '@    title "Area per lipid X Frame"'
  write(2, *) '@    xaxis  label "#Frame"'
  write(2, *) '@    yaxis  label "Area per lipid [nm\S2\N]"'

  if (rmsd)then

     inquire(file='rmsd.xvg', exist=ex)

     if (ex) then

        call execute_command_line("mv rmsd.xvg rmsd_p.xvg")
        write(*, *)
        write(*, *) 'Previous rmsd.xvg file is now backed up to rmsd_p.xvg !'

     end if

     open(3, file='rmsd.xvg', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file rmsd.xvg'
        write(*, *)
        stop

     endif

     write(3, *) '@    title "RMSD X Frame"'
     write(3, *) '@    xaxis  label "#Frame"'
     write(3, *) '@    yaxis  label "RMSD [nm]"'

  end if
  
  call cpu_time(start)

! Leitura do arquivo .pdb =================================

  i_atom = 0
  n_index = 1
  ierr = 0
  num = 1
  num2 = 1
  x_min = 1000
  y_min = 1000
  x_max = 0
  y_max = 0

  do while (ierr >= 0)

     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z
     
12   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)

     if ((atom.eq.'ATOM  '))then
        
        if(ierr > 0) then
           
           write(*, *)
           write(*, *) 'Problem reading atomic position!'
           write(*, *)
           stop
           
        endif
        
        i_atom  = i_atom + 1
        
        buff%n_atom = i_atom
        
        if (outer)then
           
           out(i_atom) = buff
           n_index = n_index + 1
           
           x_max = max(x_max, buff%x)
           x_min = min(x_min, buff%x)
           y_max = max(y_max, buff%y)
           y_min = min(y_min, buff%y)
           
        else
           
           if (i_atom == in_num(num)) then
              
              store(num) = buff
              n_index = n_index + 1
              num = num + 1
              
              x_max = max(x_max, buff%x)
              x_min = min(x_min, buff%x)
              y_max = max(y_max, buff%y)
              y_min = min(y_min, buff%y)
              
           end if
           
           if (i_atom == in_num2(num2)) then
              
              store2(num2) = buff
              n_index = n_index + 1
              num2 = num2 + 1
              
              x_max = max(x_max, buff%x)
              x_min = min(x_min, buff%x)
              y_max = max(y_max, buff%y)
              y_min = min(y_min, buff%y)
              
           end if
           
        end if
        
     else
        
        if (n_index > 1)then

           frame  = frame + 1

           eval_skip = ((mod((frame-fr_in),(n_skip+1))==0))
           
           if ((frame>fr_in-1).and.(frame<fr_end+1).and.(eval_skip)) then
              
              if (outer)then
                 
                 dx = (x_max - x_min)/bin_out
                 dy = (y_max - y_min)/bin_out
                 
                 do k=1, n_index-1
                    
                    i = nint((out(k)%x-x_min)/dx + 1)
                    
                    j = nint((out(k)%y-y_min)/dy + 1)
                    
                    if (out(k)%z>mat1(i,j)%z)then
                       
                       mat1(i,j)%z = out(k)%z
                       mat1(i,j)%x = out(k)%x
                       mat1(i,j)%y = out(k)%y
                       
                    end if
                    
                    if (out(k)%z<mat2(i,j)%z)then
                       
                       mat2(i,j)%z = out(k)%z
                       mat2(i,j)%x = out(k)%x
                       mat2(i,j)%y = out(k)%y
                       
                    end if
                    
                 end do
                 
                 n_index = 1
                 num = 1
                 num2 = 1
                 
                 do i=1,bin_out+1
                    
                    do j=1,bin_out+1
                       
                       if (mat1(i,j)%z>-1000)then
                          
                          store(num)%x = mat1(i,j)%x
                          store(num)%y = mat1(i,j)%y
                          store(num)%z = mat1(i,j)%z
                          
                          num = num + 1
                          n_index = n_index + 1
                          
                       end if
                       
                       if (mat2(i,j)%z<1000)then
                          
                          store2(num2)%x = mat2(i,j)%x
                          store2(num2)%y = mat2(i,j)%y
                          store2(num2)%z = mat2(i,j)%z
                          
                          num2 = num2 + 1
                          n_index = n_index + 1
                          
                       end if
                       
                    end do
                    
                 end do
                 
              end if! if outer   
              
              noi1 = num - 1
              noi2 = num2 - 1
              
              if (l_coarse) then
                 
                 ! Estruturação do primeiro coarse grid
                 
                 dx = (x_max - x_min)/bin_coarse
                 dy = (y_max - y_min)/bin_coarse
                 aux = 0
                 
                 !rfit relativo ao grid coarse sobre num-1 pontos
                 r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num)-1)
                 
                 al = (num-1)*100/((x_max - x_min)*(y_max - y_min))
                 al = exp(0.414754490756596*log(al)-1.29378601756983)
                 
                 do i=1, bin_coarse + 1
                    
                    do j=1, bin_coarse + 1
                       
                       aux = aux + 1
                       coarse(aux)%x = (i-1)*dx + x_min
                       coarse(aux)%y = (j-1)*dy + y_min
                       
                       s_grid = 0
                       coarse(aux)%z = 0
                       gridx = 0
                       gridy = 0
                       
                       do k=1, num - 1
                          
                          dist_x = abs(coarse(aux)%x - store(k)%x)
                          dist_y = abs(coarse(aux)%y - store(k)%y)
                          
                          dist = al*al*(dist_x**2 + dist_y**2)
                          
                          if (sqrt(dist/(al*al))<r_fit)then
                             
                             s_grid = s_grid + exp(-(dist)/pi)
                             coarse(aux)%z = coarse(aux)%z + exp(-(dist)/pi)*store(k)%z
                             gridx = gridx + exp(-(dist)/pi)*store(k)%x
                             gridy = gridy + exp(-(dist)/pi)*store(k)%y
                             
                          end if
                          
                       end do
                       
                       coarse(aux)%z = coarse(aux)%z/s_grid
                       coarse(aux)%x = gridx/s_grid
                       coarse(aux)%y = gridy/s_grid
                       
                    end do
                    
                 end do
                 
                 ! atualização do num (número de pontos do coarse grid)
                 num = aux + 1
                 
                 ! Fim da estruturação
                 ! Estruturação do segundo coarse grid
                 
                 aux = 0
                 
                 ! relativo ao grid coarse sobre num2-1 pontos
                 r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num2)-1)
                 
                 al = (num2-1)*100/((x_max - x_min)*(y_max - y_min))
                 al = exp(0.414754490756596*log(al)-1.29378601756983)
                 
                 do i=1, bin_coarse + 1
                    
                    do j=1, bin_coarse + 1
                       
                       aux = aux + 1
                       coarse2(aux)%x = (i-1)*dx + x_min
                       coarse2(aux)%y = (j-1)*dy + y_min
                       
                       s_grid = 0
                       coarse2(aux)%z = 0
                       gridx = 0
                       gridy = 0
                       
                       do k=1, num2 - 1
                          
                          dist_x = abs(coarse2(aux)%x - store2(k)%x)
                          dist_y = abs(coarse2(aux)%y - store2(k)%y)
                          
                          dist = al*al*(dist_x**2 + dist_y**2)
                          
                          if (sqrt(dist/(al*al))<r_fit)then
                             
                             coarse2(aux)%z = coarse2(aux)%z + exp(-(dist)/pi)*store2(k)%z
                             s_grid = s_grid + exp(-(dist)/pi)
                             gridx = gridx + exp(-(dist)/pi)*store2(k)%x
                             gridy = gridy + exp(-(dist)/pi)*store2(k)%y
                             
                          end if
                          
                       end do
                       
                       coarse2(aux)%z = coarse2(aux)%z/s_grid
                       coarse2(aux)%x = gridx/s_grid
                       coarse2(aux)%y = gridy/s_grid
                       
                    end do
                    
                 end do
                 
                 num2 = aux + 1
                 
              else
                 
                 coarse = store
                 coarse2 = store2
                 
              end if
              
              ! estruturação do prmeiro grid de alta resolução
              
              dx = (x_max - x_min)/n_grid
              dy = (y_max - y_min)/n_grid
              
              r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num)-1)
              
              al = (num-1)*100/((x_max - x_min)*(y_max - y_min))
              al = exp(0.414754490756596*log(al)-1.29378601756983)
              
              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    grid(i,j)%x = (i-1)*dx + x_min
                    grid(i,j)%y = (j-1)*dy + y_min
                    
                    s_grid = 0
                    grid(i,j)%z = 0
                    
                    do k=1, num - 1 
                       
                       dist_x = abs(grid(i,j)%x - coarse(k)%x)
                       dist_y = abs(grid(i,j)%y - coarse(k)%y)
                       
                       dist = al*al*(dist_x**2 + dist_y**2)
                       
                       if (sqrt(dist/(al*al))<r_fit)then
                          
                          s_grid = s_grid + exp(-(dist)/pi)
                          grid(i,j)%z = grid(i,j)%z + exp(-(dist)/pi)*coarse(k)%z
                          
                       end if
                       
                    end do
                    
                    grid(i,j)%z = grid(i,j)%z/s_grid
                    
                 end do
                 
              end do
              
              
              ! estruturação do segundo grid de alta resolução
              
              r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(num)-1)
              
              al = (num-1)*100/((x_max - x_min)*(y_max - y_min))
              al = exp(0.414754490756596*log(al)-1.29378601756983)
              
              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    grid2(i,j)%x = (i-1)*dx + x_min
                    grid2(i,j)%y = (j-1)*dy + y_min
                    
                    s_grid = 0
                    grid2(i,j)%z = 0
                    
                    do k=1, num2 - 1
                       
                       dist_x = abs(grid2(i,j)%x - coarse2(k)%x)
                       dist_y = abs(grid2(i,j)%y - coarse2(k)%y)
                       
                       dist = al*al*(dist_x**2 + dist_y**2)
                       
                       if (sqrt(dist/(al*al))<r_fit)then
                          
                          s_grid = s_grid + exp(-(dist)/pi)
                          grid2(i,j)%z = grid2(i,j)%z + exp(-(dist)/pi)*coarse2(k)%z
                          
                       end if
                       
                    end do
                    
                    grid2(i,j)%z = grid2(i,j)%z/s_grid
                    
                 end do
                 
              end do
              
              ! Fim da estruturação
              
              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    grid3(i,j)%z = (grid(i,j)%z + grid2(i,j)%z)/2
                    grid3(i,j)%x = grid(i,j)%x
                    grid3(i,j)%y = grid(i,j)%y
                    
                 end do
                 
              end do
              
              if (rmsd)then
                 
                 noir = 0
                 
                 do i=1, noi1
                    
                    a = nint((store(i)%x-x_min)/dx) + 1
                    
                    b = nint((store(i)%y-y_min)/dy) + 1
                    
                    dist_z = store(i)%z-grid(a,b)%z
                    
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 do i=1, noi2
                    
                    a = nint((store2(i)%x-x_min)/dx) + 1
                    
                    b = nint((store2(i)%y-y_min)/dy) + 1
                    
                    dist_z = store2(i)%z-grid2(a,b)%z
                    
                    noir = noir + dist_z*dist_z/(noi1+noi2)
                    
                 end do
                 
                 noir = sqrt(noir)
                 
                 write(3, *) noir/10
                 
              end if
              
              ! Calculando a área pela fórmula de Heron
              
              s_area = 0
              
              do i=2, n_grid+1
                 
                 do j=2, n_grid+1
                    
                    dist_x = abs(grid3(i-1,j-1)%x - grid3(i-1,j)%x)
                    dist_y = abs(grid3(i-1,j-1)%y - grid3(i-1,j)%y)
                    dist_z = abs(grid3(i-1,j-1)%z - grid3(i-1,j)%z)
                    
                    la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i-1,j-1)%x - grid3(i,j-1)%x)
                    dist_y = abs(grid3(i-1,j-1)%y - grid3(i,j-1)%y)
                    dist_z = abs(grid3(i-1,j-1)%z - grid3(i,j-1)%z)
                    
                    lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i-1,j)%x - grid3(i,j-1)%x)
                    dist_y = abs(grid3(i-1,j)%y - grid3(i,j-1)%y)
                    dist_z = abs(grid3(i-1,j)%z - grid3(i,j-1)%z)
                    
                    lc = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    s_soma = (la + lb + lc)/2
                    
                    s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
                    
                    ! segundo triângulo
                    
                    dist_x = abs(grid3(i-1,j)%x - grid3(i,j)%x)
                    dist_y = abs(grid3(i-1,j)%y - grid3(i,j)%y)
                    dist_z = abs(grid3(i-1,j)%z - grid3(i,j)%z)
                    
                    la = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    dist_x = abs(grid3(i,j-1)%x - grid3(i,j)%x)
                    dist_y = abs(grid3(i,j-1)%y - grid3(i,j)%y)
                    dist_z = abs(grid3(i,j-1)%z - grid3(i,j)%z)
                    
                    lb = sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                    
                    s_soma = (la + lb + lc)/2
                    
                    s_area = s_area + sqrt(s_soma*(s_soma-la)*(s_soma-lb)*(s_soma-lc))
                    
                 end do
                 
              end do
              
              ! Fim do cálculo
              
              write(2, *) frame, s_area/100/n_lipid
              
           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
           
           n_index = 1
           i_atom = 0
           num = 1
           num2 = 1
           x_min = 1000
           y_min = 1000
           x_max = 0
           y_max = 0

           !====garante que fr_end sempre seja maior que frame ===
           !====caso essa variável não tenha sido fixada==========
           if (.not.end)then
              
              fr_end = frame + 1
              
           end if
           !======================================================
           
           if (outer)then
              
              do i=1, bin_out + 1
                 
                 do j=1, bin_out + 1
                    
                    mat1(i,j)%z = -1000
                    mat2(i,j)%z = 1000
                    
                 end do
                 
              end do
              
           end if
           
        end if
        
     end if
     
  end do
  
  call cpu_time(finish)
  
  ! Fim da leitura do arquivo .pdb

  close(1)
  close(2)

  if (rmsd)then
     
     close(3)

  end if

  if (p_grid)then

     inquire(file='grid1.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid1.pdb grid1_p.pdb")
        write(*, *)
        write(*, *) 'Previous grid1.pdb file is now backed up to grid1_p.pdb !'

     end if

     open(3, file='grid1.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid1.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid(i,j)%x, &
                grid(i,j)%y, grid(i,j)%z

        end do

     end do

     close(3)
     !first grid created

     inquire(file='grid2.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid2.pdb grid2_p.pdb")
        write(*, *)
        write(*, *) 'Previous grid2.pdb file is now backed up to grid2_p.pdb !'

     end if

     open(3, file='grid2.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid2.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid2(i,j)%x, &
                grid2(i,j)%y, grid2(i,j)%z

        end do

     end do

     close(3)
     !second grid created

     inquire(file='grid3.pdb', exist=ex)

     if (ex) then

        call execute_command_line("mv grid3.pdb grid3_p.pdb")
        write(*, *)
        write(*, *) 'Previous grid3.pdb file is now backed up to grid3_p.pdb !'

     end if

     open(3, file='grid3.pdb', status='new', iostat=ierr)

     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file grid3.pdb'
        write(*, *)
        stop

     end if

     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'

     do i=1, n_grid+1

        do j=1, n_grid+1

           write(3, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                '  DOT',' ', 1,'    ', grid3(i,j)%x, &
                grid3(i,j)%y, grid3(i,j)%z

        end do

     end do

     close(3)
     !third grid created
     
  end if

  write(*, *)
  write(*, *) "Finished"

  
  
  hour = (finish-start)/3600
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)
  
  if (outer)then

     deallocate(mat1)
     deallocate(mat2)

  end if
  
end program ar_grid
