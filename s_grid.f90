module grd

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


end module grd

program grid_project

  use grd

  logical :: ex, bin, outer, down, rmsd, l_coarse
  logical :: begin, end, skip, eval_skip
  
  integer :: i, j, k, ierr, n_grid, frame, aux, num, noi1, noi2 
  integer :: n_index, n_atom, i_atom, bin_out, bin_coarse, a, b
  integer :: fr_in, fr_end, n_skip
  integer, dimension(50000) :: in_num

  real :: dx, dy, dist, s_grid, r_fit, hour, minu, sec, start, finish
  real :: x_max, x_min, y_max, y_min, aux2, gridx, gridy, noir, al
  real :: z_max, z_min

  character(len=30) :: coord,  index, atom
  character(len=30), dimension(20) :: get

  type(vet1) :: buff
  type(vet1), dimension(50000) :: store, coarse
  type(vet1), dimension(500000) :: out
  type(vet2), dimension(:,:), allocatable :: mat
  type(vet2), dimension(1001,1001) :: grid

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
  write(*, *) "                       ** s_grid **"
  write(*, *) ""
  write(*, *) "             Santos, D. E. S.; Soares, T. A."
  write(*, *) ""
  write(*, *) "Please cite SuAVE: A Tool for Analyzing Curvature-Dependent" 
  write(*, *) "Properties in Chemical Interfaces. Denys E. S. Santos," 
  write(*, *) "Frederico J. S. Pontes, Roberto D. Lins, Kaline Coutinho," 
  write(*, *) "Thereza A. Soares. J. Chem. Inf. Model. 2019."
  write(*, *) ""
  write(*, *) ""
  write(*, *) "s_grid builds a grid per frame throughout a trajectory file." 
  write(*, *) "Its output is useful to verify how accurate is the fitting of" 
  write(*, *) "the calculated grid on the chemical surface.  "
  write(*, *) ""
  write(*, *) "Usage: s_grid -in file.pdb -ind file.ndx"
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
  write(*, *) ""
  Write(*, *) "-down used to generate a PDB file with the last "
  Write(*, *) "external lowest fitting grid. This option is only used when" 
  Write(*, *) "followed by option –outer."  

!
! pegando os arquivos de entrada
!

  outer = .false.
  bin = .false.
  down = .false.
  coord = 'miss'
  index = 'miss'
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

     if (get(i)=='-ind')then

        index = get(i+1)

     end if

     if (get(i)=='-bin')then

        bin= .true.
        
     end if

     if (get(i)=='-outer')then

        outer= .true.

     end if

     if (get(i)=='-down')then

        down= .true.

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

     if (index=='miss')then

        write(*, *)
        write(*, *)'Index file is missing'
        write(*, *)
        stop

     end if

  end if
  
  !==Lendo os arquivos de index=========================
  if (.not.outer)then
     
     open(2, file=index, status='old', iostat=ierr)
     
     if(ierr /= 0) then

        write(*, *)
        write(*, *) 'Unable to open file', index
        write(*, *)
        stop
        
     end if

     ! Inicio da leitura do index ============================
     
     n_index = 1
     
     do while (ierr==0)
        
        read(2, *, iostat=ierr) aux
        
        if (ierr == 0)then
           
           in_num(n_index) = aux
           n_index = n_index + 1
           
        end if
        
     end do
     
     n_atom  = n_index - 1
     
     if(ierr>0) then

        write(*, *)
        write(*, *) 'Problem by reading ', index
        write(*, *)
        stop
        
     end if
     
     do i=n_index, 50000
        
        in_num(i)=-1
        
     end do

     write(*, *)
     write(*, *) index, " input file is OK"
     
     close(2)

  end if
  
  
2 format(a10)

  if ((bin).or.(outer)) then

     write(*, *)
     write(*, 2, advance='no') '    bin = '
     read(*, *, iostat=ierr) n_grid
     
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
    
  !=======================================================
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
  
  !==================================================
  
  !==================================================
  ! definição da matriz de pontos externos
  
  if (outer)then

     allocate (mat(bin_out+1, bin_out+1), stat = ierr)
     if (ierr /= 0) stop 'Not enough memory to initialize mat matrix'
     
     do i=1, bin_out + 1
        
        do j=1, bin_out + 1

           if (down)then
              
              mat(i,j)%z = 1000

           else

              mat(i,j)%z = -1000

           end if

        end do
        
     end do

  end if

  !==================================================

  write(*, *)
  write(*, *) 'Processing ..........' !==============================iniciando processamento

  open(1, file=coord, status='old', iostat=ierr)

  if(ierr /=0) then

     write(*, *)
     write(*, *) 'Unable to open file ', coord
     write(*, *)
     stop

  endif

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

  !==========================================================
  inquire(file='grid.pdb', exist=ex)

  if (ex) then

     call execute_command_line("mv grid.pdb grid_p.pdb")
     write(*, *)
     write(*, *) 'Previous grid.pdb file is now backed up to grid_p.pdb !'

  end if

  open(4, file='grid.pdb', status='new', iostat=ierr)

  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file grid.pdb'
     write(*, *)
     stop

  end if

  call cpu_time(start)
  
  ! Leitura do arquivo .pdb
  
  i_atom = 0
  frame = 0
  ierr = 0
  n_index = 1
  x_min = 1000
  y_min = 1000
  z_min = 1000
  x_max = 0
  y_max = 0
  z_max = 0

  do while (ierr>=0)
  
     read(1, 12, iostat=ierr) atom, buff%n_atom, buff%atom, &
          buff%resid, buff%ident, buff%n_resid, buff%code, buff%x, &
          buff%y, buff%z

12   format(a6, i5.1, a5, a5, a1, i4.1, a4, 3f8.3)
     
     if (atom.eq.'ATOM  ') then

        if(ierr > 0) then

           write(*, *)
           write(*, *) 'Problem reading atomic position!'
           write(*, *)
           stop

        endif

        i_atom = i_atom + 1

        buff%n_atom = i_atom
        
        if (outer)then

           out(i_atom) = buff
           n_index = n_index + 1

           x_max = max(x_max, buff%x)
           x_min = min(x_min, buff%x)
           y_max = max(y_max, buff%y)
           y_min = min(y_min, buff%y)
           z_max = max(z_max, buff%z)
           z_min = min(z_min, buff%z)

        else
           
           if (i_atom == in_num(n_index)) then
              
              store(n_index) = buff
              n_index = n_index + 1

              x_max = max(x_max, buff%x)
              x_min = min(x_min, buff%x)
              y_max = max(y_max, buff%y)
              y_min = min(y_min, buff%y)
              z_max = max(z_max, buff%z)
              z_min = min(z_min, buff%z)

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
                    
                    if (down) then
                       
                       if (out(k)%z<mat(i,j)%z)then
                          
                          mat(i,j)%z = out(k)%z
                          mat(i,j)%x = out(k)%x
                          mat(i,j)%y = out(k)%y
                          
                       end if
                       
                    else 
                       
                       if (out(k)%z>mat(i,j)%z)then
                          
                          mat(i,j)%z = out(k)%z
                          mat(i,j)%x = out(k)%x
                          mat(i,j)%y = out(k)%y
                          
                       end if
                       
                    end if
                    
                 end do
                 
                 num = n_index
                 n_index = 1
                 
                 do i=1,bin_out+1
                    
                    do j=1,bin_out+1
                       
                       if (down)then
                          
                          if (mat(i,j)%z<1000)then
                             
                             store(n_index)%x = mat(i,j)%x
                             store(n_index)%y = mat(i,j)%y
                             store(n_index)%z = mat(i,j)%z
                             
                             n_index = n_index + 1
                             
                          end if
                          
                       else
                          
                          if (mat(i,j)%z>-1000)then
                             
                             store(n_index)%x = mat(i,j)%x
                             store(n_index)%y = mat(i,j)%y
                             store(n_index)%z = mat(i,j)%z
                             
                             n_index = n_index + 1
                             
                          end if
                          
                       end if
                       
                    end do
                    
                 end do
                 
              end if! if outer
              
              noi1 = n_index - 1
              
              if (l_coarse) then
                 
                 !
                 ! Estruturação do grid
                 
                 aux = 0
                 dx = (x_max - x_min)/bin_coarse
                 dy = (y_max - y_min)/bin_coarse
                 
                 ! refinamento do r_fit e do alfa
                 r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(n_index)-1)
                 
                 al = (n_index-1)*100/((x_max - x_min)*(y_max - y_min))
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
                       
                       do k=1, n_index - 1
                          
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
                 
                 n_index = aux + 1
                 
              else
                 
                 coarse = store
                 
              end if
              
              ! estruturação do prmeiro grid de alta resolução
              dx = (x_max - x_min)/n_grid
              dy = (y_max - y_min)/n_grid
              
              ! refinamento do r_fit e do alfa
              r_fit = 3*sqrt((x_max - x_min)**2 + (y_max - y_min)**2)/sqrt(real(n_index)-1)
              
              al = (n_index-1)*100/((x_max - x_min)*(y_max - y_min))
              al = exp(0.414754490756596*log(al)-1.29378601756983)
              
              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    grid(i,j)%x = (i-1)*dx + x_min
                    grid(i,j)%y = (j-1)*dy + y_min
                    
                    s_grid = 0
                    grid(i,j)%z = 0
                    
                    do k=1, n_index - 1
                       
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
              
              if (rmsd)then
                 
                 noir = 0
                 
                 do i=1, noi1
                    
                    a = nint((store(i)%x-x_min)/dx) + 1
                    
                    b = nint((store(i)%y-y_min)/dy) + 1
                    
                    dist_z = store(i)%z-grid(a,b)%z
                    
                    noir = noir + dist_z*dist_z/noi1
                    
                 end do
                 
                 noir = sqrt(noir)
                 
                 write(3, *) noir/10
                 
              end if
              
              ! escreve cada frame da trajetoria do grid

3             format(a8, f7.2, a2, f7.2, a2, f7.2, a35)
              write(4, '(a26)') 'REMARK Generated by s_grid'
              write(4, '(a34)') 'TITLE     Rectangular Regular GRID'
              write(4, '(a22)') 'REMARK trajectory file'
              write(4, 3) 'CRYST1  ', x_max-x_min, '  ', y_max-y_min, '  ', z_max-z_min, '  90.00  90.00  90.00 P  1       1 '
              write(4, '(a6, i4)') 'MODEL ', frame

              do i=1, n_grid+1
                 
                 do j=1, n_grid+1
                    
                    write(4, 12, iostat=ierr) 'ATOM  ', (i-1)*(n_grid+1) + j, '  DOT', &
                         '  DOT',' ', 1,'    ', grid(i,j)%x, &
                         grid(i,j)%y, grid(i,j)%z
                    
                 end do
                 
              end do
              
              write(4, '(a3)') 'TER'
              write(4, '(a6)') 'ENDMDL'
              
           !fim da escrita da trajetoria

           end if !======((frame<fr_in-1).and.(frame>fr_end+1))
              
           n_atom = n_index-1
           n_index = 1
           i_atom = 0
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

                    if (down)then

                       mat(i,j)%z = 1000

                    else

                       mat(i,j)%z = -1000

                    end if
                    
                 end do
                 
              end do

           end if
           
        end if
        
     end if 
     
  end do
  
  ! Fim da estruturação

  close(1)
  close(4)

  if (rmsd)then
     
     close(3)
     
  end if
  
  call cpu_time(finish)

  
  ! desenho dos pontos usados no ajuste
  
  inquire(file='adjust.pdb', exist=ex)
  
  if (ex) then
     
     call execute_command_line("mv adjust.pdb adjust_p.pdb")
     write(*, *)
     write(*, *) 'Previous adjust.pdb file is now backed up to adjust_p.pdb !'
     
  end if
  
  open(3, file='adjust.pdb', status='new', iostat=ierr)
  
  if(ierr /= 0) then

     write(*, *)
     write(*, *) 'Unable to open file ajuste.pdb'
     write(*, *)
     stop

  end if

  write(3, *) 'TITLE     Rectangular Regular GRID'
  write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
  write(3, *) 'MODEL        1'

  do i=1, noi1
     
     write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
          '  DOT',' ', 1,'    ', store(i)%x, &
          store(i)%y, store(i)%z

  end do

  close(3)
  
  if (l_coarse) then
     
     inquire(file='coarse.pdb', exist=ex)
     
     if (ex) then
        
        call execute_command_line("mv coarse.pdb coarse_p.pdb")
        write(*, *)
        write(*, *) 'Previous coarse.pdb file is now backed up to coarse_p.pdb !'
        
     end if
     
     open(3, file='coarse.pdb', status='new', iostat=ierr)
     
     if(ierr /= 0) then
        
        write(*, *)
        write(*, *) 'Unable to open file coarse.pdb'
        write(*, *)
        stop
        
     end if
     
     write(3, *) 'TITLE     Rectangular Regular GRID'
     write(3, *) 'REMARK    THIS IS A SIMULATION BOX'
     write(3, *) 'MODEL        1'
     
     do i=1, (bin_coarse + 1)*(bin_coarse + 1)
        
        write(3, 12, iostat=ierr) 'ATOM  ', i, '  DOT', &
             '  DOT',' ', 1,'    ', coarse(i)%x, &
             coarse(i)%y, coarse(i)%z
        
     end do

  end if
  
  write(*, *)
  write(*, *)'Finished'
  
! Fim. Arquivo escrito

  close(3)

  hour = (finish-start)/3600
  minu = (hour - int(hour))*60
  sec = (minu-int(minu))*60

  write(*, *)
  write(*, 1) " Processing time : ", int(hour), " h", int(minu), " min", int(sec), " sec"
1 format(a19, i6.3, a2, i6.3, a4, i6.3, a4)
  write(*, *)
  
  if (outer)then
     
     deallocate(mat)

  end if
     
end program grid_project
