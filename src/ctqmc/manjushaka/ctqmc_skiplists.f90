!>>> this module define skip lists 
  module m_skiplists
     use constants
     use control
     use context
     use m_sector

     implicit none

!>>> data structure of skip lists 
     type :: t_skiplists
         integer :: list_level(2)
         integer :: num_node(2)
         type(t_node), pointer :: head
         type(t_node), pointer :: tail
     end type t_skiplists

!>>> data structure of a node
     type :: t_node
         integer :: node_level(2)
         integer :: indx(2)
         type(t_pnode), dimension(:,:), pointer     :: forward => null()
         logical, dimension(:,:,:), pointer         :: lsave   => null() 
         logical, dimension(:,:), pointer           :: lcopy => null()
         type(t_subprod), dimension(:,:), pointer   :: subprod => null()
     end type t_node

!>>> data structure of a pointer to node
     type :: t_pnode
         type(t_node), pointer :: p
     end type t_pnode

!>>> data structure of a subprod
     type :: t_subprod
         real(dp), dimension(:,:,:), pointer :: submat => null()
     end type t_subprod

!>>> a very large integer
     integer, private, parameter :: max_int = 100000000

!>>> global variables, we use two skip lists
     type(t_skiplists), public, save, pointer :: skip_lists
     type(t_node), public, save, pointer :: node_a
     type(t_node), public, save, pointer :: node_b

     contains

!>>> create a new skiplists
     subroutine new_skiplists(skiplists)
        implicit none

        type(t_skiplists), pointer :: skiplists
  
        type(t_node), pointer :: head, tail
        integer :: i

! allocate memory for skiplists
        allocate(skiplists)

! add a head node with maximum level, and a tail node with minimum level
        call new_node(mlevl, head)
        call new_node(1,     tail)

! set head node
        head%indx = 0
        do i=1, mlevl
            head%forward(i,1)%p => tail
            head%forward(i,2)%p => tail
        enddo
        head%lsave = .false. 
        head%lcopy = .false.

! set tail node
        tail%indx = max_int
        tail%forward(1,1)%p => null()
        tail%forward(1,2)%p => null()
 
! set skiplists
        skiplists%list_level = 2
        skiplists%num_node = 2
        skiplists%head => head
        skiplists%tail => tail

        return
     end subroutine new_skiplists

!>>> destroy a skiplists
     subroutine destroy_skiplists(skiplists)
        implicit none

        type(t_skiplists), pointer :: skiplists
 
        type(t_pnode) :: curr, next

        if (associated(skiplists)) then
            curr%p => skiplists%head
            do while(associated(curr%p))
                next%p => curr%p%forward(1,2)%p
                call destroy_node(curr%p)
                curr%p => next%p
            enddo
            deallocate(skiplists)
        endif
       
        return
     end subroutine destroy_skiplists

!>>> crate a new node
     subroutine new_node(level, node)
        implicit none
        
        integer, intent(in) :: level 
        type(t_node), pointer :: node
        integer :: i,j

        allocate(node)

        node%node_level = level
        if (level >= 1) then
            allocate( node%forward(level, 2) )
            allocate( node%lsave(nsectors, level, 2) )
            allocate( node%lcopy(nsectors, level) )
            allocate( node%subprod(nsectors, level) ) 
            do i=1, level
                do j=1, nsectors
                    call new_subprod_item(node%subprod(j,i))
                enddo
            enddo
            node%lsave = .false.
            node%lcopy = .false.
        endif
     
        return
     end subroutine new_node

!>>> destroy a node
     subroutine destroy_node(node)
        implicit none

        type(t_node), pointer :: node

        integer :: i,j

        if (associated(node)) then
            if (associated(node%forward))   deallocate(node%forward) 
            if (associated(node%lsave))     deallocate(node%lsave)
            if (associated(node%lcopy))     deallocate(node%lcopy)
            if (associated(node%subprod)) then
                do i=1, node%node_level(1)
                    do j=1, nsectors
                        call destroy_subprod_item(node%subprod(j,i))
                        call destroy_subprod_item(node%subprod(j,i))
                    enddo
                enddo
                deallocate(node%subprod)
            endif

            deallocate(node) 
        endif
   
        return
     end subroutine destroy_node

!>>> create a new subprod
     subroutine new_subprod_item(subprod)
        implicit none

        type(t_subprod), intent(inout) :: subprod

        if ( .not. associated(subprod%submat) ) then
            allocate( subprod%submat(max_dim_sect, max_dim_sect, 2) )
        endif 

        return
     end subroutine new_subprod_item

!>>> destroy a subprod
     subroutine destroy_subprod_item(subprod)
        implicit none

        type(t_subprod), intent(inout) :: subprod

        if (associated(subprod%submat)) then
            deallocate(subprod%submat)
        endif

        return
     end subroutine destroy_subprod_item

!>>> insert a node 
     subroutine insert_node(skiplists, level, indx, node)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(inout) :: level
        integer, intent(in) :: indx
        type(t_node), pointer :: node

        type(t_pnode) :: update(mlevl), x
        integer :: i,j

        x%p => skiplists%head
        do i=skiplists%list_level(1), 1, -1
            do while (x%p%forward(i,1)%p%indx(1) <= indx)
                x%p => x%p%forward(i,1)%p
            enddo 
            x%p%lsave(:,i,1) = .false.
            update(i)%p => x%p
        enddo 

        if (level > skiplists%list_level(1) - 1) then
            level = skiplists%list_level(1)
            update(level)%p => skiplists%head
            skiplists%head%lsave(:,level,1) = .false.

            skiplists%list_level(1) = skiplists%list_level(1) + 1
            skiplists%head%forward(skiplists%list_level(1),1)%p => skiplists%tail
            skiplists%head%lsave(:, skiplists%list_level(1),1) = .false.
        endif 

        call new_node(level, node)
        node%indx(1) = indx+1
        
        do i=level, 1, -1
            node%lsave(:,i,1) = .false.
            x = update(i)
            node%forward(i,1)%p => x%p%forward(i,1)%p
            x%p%forward(i,1)%p => node
        enddo

        x%p => node%forward(1,1)%p
        do while( .not. associated(x%p, skiplists%tail) )
            x%p%indx(1) = x%p%indx(1) + 1
            x%p => x%p%forward(1,1)%p
        enddo

        skiplists%num_node(1) = skiplists%num_node(1) + 1
  
        return
     end subroutine insert_node

!>>> remove a node
     subroutine remove_node(skiplists, indx, node)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: indx
        type(t_node), pointer :: node

        type(t_pnode) :: update(mlevl), x
        integer :: level
        integer :: i, j

        x%p => skiplists%head
        do i=skiplists%list_level(1), 1, -1
            do while (x%p%forward(i,1)%p%indx(1) < indx)
                x%p => x%p%forward(i,1)%p
            enddo 
            x%p%lsave(:,i,1) = .false.
            update(i)%p => x%p
        enddo 

        node => x%p%forward(1,1)%p

        x%p => node%forward(1,1)%p
        do while(.not. associated(x%p, skiplists%tail))
            x%p%indx(1) = x%p%indx(1) - 1
            x%p => x%p%forward(1,1)%p
        enddo

        do i=1, node%node_level(1)
            update(i)%p%forward(i,1)%p => node%forward(i,1)%p 
        enddo

        do while (skiplists%list_level(1)>2) 
            level = skiplists%list_level(1) - 1
            x%p => skiplists%head%forward(level, 1)%p
            if ( associated( x%p, skiplists%tail ) ) then
                skiplists%list_level(1) = level
            else
                EXIT 
            endif
        enddo

        skiplists%num_node(1) = skiplists%num_node(1) - 1

        return
     end subroutine remove_node

!>>> change a node
     subroutine change_one_node(skiplists, indx)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: indx

        type(t_pnode) :: x
        integer :: i,j

        x%p => skiplists%head
        do i=skiplists%list_level(1), 1, -1
            do while (x%p%forward(i,1)%p%indx(1) <= indx)
                x%p => x%p%forward(i,1)%p
            enddo 
            x%p%lsave(:,i,1) = .false.
        enddo 

        return
     end subroutine change_one_node

!>>> change all node
     subroutine change_all_node(skiplists)
        implicit none

        type(t_skiplists), pointer :: skiplists

        type(t_pnode) :: x
        integer :: i,j
 
        do i=1, skiplists%list_level(1)
            x%p => skiplists%head
            do while(.not. associated(x%p, skiplists%tail))
                x%p%lsave(:,i,1) = .false.
                x%p => x%p%forward(i,1)%p
            enddo
        enddo

        return
     end subroutine change_all_node 

!>>> deep copy a skip lists
     subroutine deep_copy_skiplists(skiplists, copy_mode)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: copy_mode

        type(t_pnode) :: x
        integer :: i,j

        if (copy_mode == 1) then
            skiplists%list_level(2) = skiplists%list_level(1)
            skiplists%num_node(2)   = skiplists%num_node(1)
            x%p => skiplists%head
            do while(.not. associated(x%p, skiplists%tail))
                x%p%node_level(2) = x%p%node_level(1)
                x%p%indx(2) = x%p%indx(1)
                x%p%forward(:,2) = x%p%forward(:,1)
                x%p%lsave(:,:,2) = x%p%lsave(:,:,1)
                do i=1, x%p%node_level(1)
                    do j=1, nsectors
                        if (x%p%lcopy(j,i)) then
                            x%p%subprod(j,i)%submat(:,:,2) = x%p%subprod(j,i)%submat(:,:,1)
                        endif                 
                    enddo
                enddo
                x%p => x%p%forward(1,1)%p
            enddo
        else
            skiplists%list_level(1) = skiplists%list_level(2)
            skiplists%num_node(1)   = skiplists%num_node(2)
            x%p => skiplists%head
            do while(.not. associated(x%p, skiplists%tail))
                x%p%node_level(1) = x%p%node_level(2)
                x%p%indx(1) = x%p%indx(2)
                x%p%forward(:,1) = x%p%forward(:,2)
                x%p%lsave(:,:,1) = x%p%lsave(:,:,2)
                x%p%lcopy(:,:) = .false.
                x%p => x%p%forward(1,1)%p
            enddo
        endif

        return
     end subroutine deep_copy_skiplists

     subroutine trial_update_skiplists(skiplists, csize, imove, ia, ib)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: csize
        integer, intent(in) :: imove
        integer, intent(in) :: ia
        integer, intent(in) :: ib

        integer :: level

! first, copy 2 to 1
        call deep_copy_skiplists(skiplists, 2)

        if (imove == 1) then
            call random_level(level)
            call insert_node(skiplists, level, ia-1, node_a) 
            if (ia < csize-1) then
                call change_one_node(skiplists,ia+1)
            endif
            call random_level(level)
            call insert_node(skiplists, level, ib-1, node_b) 
            if (ib < csize) then
                call change_one_node(skiplists,ib+1)
            endif
! remove two operators
        elseif (imove == 2) then
            call remove_node(skiplists, ia, node_a) 
            if (ia < csize+2) then
                call change_one_node(skiplists,ia)
            endif
            call remove_node(skiplists, ib, node_b) 
            if (ib < csize+1) then
                call change_one_node(skiplists,ib)
            endif
! remove one, then insert one
        elseif (imove == 3 .or. imove == 4) then
            call remove_node(skiplists, ia, node_a) 
            if (ia < csize) then
                call change_one_node(skiplists,ia)
            endif
            call random_level(level)
            call insert_node(skiplists, level, ib-1, node_b) 
            if (ib < csize) then
                call change_one_node(skiplists,ib+1)
            endif
! change all the node
        elseif (imove == 5) then
            call change_all_node(skiplists)
        endif

        return
     end subroutine trial_update_skiplists

     subroutine real_update_skiplists(skiplists, imove, accept)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: imove
        logical, intent(in) :: accept

        if (accept) then
            if (imove == 2) then
                call destroy_node(node_a)
                call destroy_node(node_b)
            elseif(imove == 3 .or. imove == 4) then
                call destroy_node(node_a)
            endif
            call deep_copy_skiplists(skiplists, 1)
        else
            if (imove == 1) then
                call destroy_node(node_a)
                call destroy_node(node_b)
            elseif (imove == 3 .or. imove == 4) then
                call destroy_node(node_b)
            endif
        endif
        
        return
     end subroutine real_update_skiplists

     subroutine ctqmc_sector_ztrace(skiplists, csize, string, index_t_loc, expt_t_loc, trace)
        implicit none
     
        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: csize
        integer, intent(in) :: string(csize+1)
        integer, intent(in) :: index_t_loc(mkink)
        real(dp), intent(in) :: expt_t_loc(ncfgs)
        real(dp), intent(out) :: trace
        
        real(dp) :: right_mat(max_dim_sect, max_dim_sect)
        type(t_node), pointer :: head
        type(t_node), pointer :: tail
        integer :: isect
        integer :: level
        integer :: indx
        integer :: dim1
        integer :: i,j
     
        level = skiplists%list_level(1)
        head => skiplists%head
        tail => skiplists%tail
     
        isect = string(1)
        dim1 = sectors(isect)%ndim

! if we don't save its result, allocate memory and calculate it
        !if (.not. head%lsave(isect, level, 1)) then
        !    call cat_subprod(skiplists, level, head, tail, csize, string, index_t_loc, right_mat, head%subprod(isect,level)%submat(:,:,1))
        !    head%lsave(isect, level, 1) = .true.
        !    head%lcopy(isect, level) = .true.
        !    right_mat = head%subprod(isect, level)%submat(:,:,1)
        !else
        !    right_mat = head%subprod(isect, level)%submat(:,:,2)
        !endif

        call cat_subprod2(skiplists, csize, string, index_t_loc, right_mat)

! special treatment of the last time-evolution operator
        indx = sectors(isect)%istart
        do i=1, dim1
            do j=1, dim1
                right_mat(i,j) = right_mat(i,j) * expt_t_loc(indx+i-1)
            enddo
        enddo
        num_prod = num_prod + one

! store final product
        sectors(isect)%final_product(:,:,1) = right_mat(1:dim1, 1:dim1)

! calculate the trace
        trace  = zero
        do j=1, dim1
            trace = trace + right_mat(j,j)
        enddo

        return
     end subroutine ctqmc_sector_ztrace
     
     recursive &
     subroutine cat_subprod(skiplists, level, head, tail, csize, string, index_t_loc, tmp_mat, total_prod)
        implicit none
     
        type(t_skiplists), pointer :: skiplists
        integer, intent(in)   :: level
        type(t_node), pointer :: head
        type(t_node), pointer :: tail
        integer, intent(in)   :: csize
        integer, intent(in)   :: string(csize+1)
        integer, intent(in)   :: index_t_loc(mkink)
        real(dp), intent(inout) :: total_prod(max_dim_sect, max_dim_sect)
        real(dp), intent(out) :: tmp_mat(max_dim_sect, max_dim_sect)

        type(t_pnode) :: x1, x2
        integer :: this_level
        integer :: counter
        integer :: isect, jsect
        integer :: dim1, dim2, dim3
        integer :: istart
        integer :: vt, vf
        integer :: i,j
     
        this_level = level - 1
        total_prod = zero

        if (head%indx(1) == 0) then
            dim1 = sectors(string(1))%ndim
        else
            dim1 = sectors(string(head%indx(1)))%ndim
        endif
 
        counter = 0 
        x1%p => head
        x2%p => head%forward(this_level,1)%p

        if (this_level == 1) then
            do while(.not. associated(x1%p, tail))
                counter = counter + 1
                if (x1%p%indx(1) == 0) then
                    isect = string(1) 
                    jsect = string(1) 
                else
                    isect = string(x1%p%indx(1))
                    jsect = string(x1%p%indx(1)+1)
                endif

                dim2 = sectors(isect)%ndim
                dim3 = sectors(jsect)%ndim

                if (.not. x1%p%lsave(isect,1,1)) then
                    x1%p%subprod(isect,1)%submat(:,:,1) = zero
                    if (x1%p%indx(1) == 0) then
                        do i=1, dim2
                            x1%p%subprod(isect,1)%submat(i,i,1) = one 
                        enddo
                    else
                        istart = sectors(isect)%istart
                        vt = type_v( index_t_loc(x1%p%indx(1)) )
                        vf = flvr_v( index_t_loc(x1%p%indx(1)) ) 
                        do i=1, dim2
                            do j=1, dim3
                                x1%p%subprod(isect,1)%submat(j,i,1) = sectors(isect)%myfmat(vf,vt)%item(j,i) * &
                                expt_v( istart+i-1, index_t_loc(x1%p%indx(1)) )
                            enddo 
                        enddo
                        num_prod = num_prod + 1
                    endif
                    x1%p%lsave(isect,1,1) = .true.
                    x1%p%lcopy(isect,1) = .true.

                    if (counter == 1) then
                        total_prod = x1%p%subprod(isect,1)%submat(:,:,1)
                    else
                        call dgemm( 'N','N', dim3, dim1, dim2, one,                    &
                                    x1%p%subprod(isect,1)%submat(:,:,1), max_dim_sect, &
                                    total_prod,                          max_dim_sect, &
                                    zero, tmp_mat,                       max_dim_sect   )
                        total_prod = tmp_mat
                        num_prod = num_prod + 1
                    endif

                else
                    if (counter == 1) then
                        total_prod = x1%p%subprod(isect,1)%submat(:,:,2)
                    else
                        call dgemm( 'N','N', dim3, dim1, dim2, one,                  &
                                    x1%p%subprod(isect,1)%submat(:,:,2), max_dim_sect, &
                                    total_prod,                        max_dim_sect, &
                                    zero, tmp_mat,                     max_dim_sect   )
                        total_prod = tmp_mat
                        num_prod = num_prod + 1
                    endif
                endif

                x1%p => x1%p%forward(1,1)%p
            enddo

        else
            do while(.not. associated(x1%p, tail) )
                counter =  counter + 1
                if (x1%p%indx(1) == 0) then
                    isect = string(1) 
                    if (associated(x2%p, skiplists%tail)) then
                        jsect = string(1)
                    else
                        jsect = string(x2%p%indx(1)) 
                    endif
                else
                    isect = string(x1%p%indx(1))
                    if (associated(x2%p, skiplists%tail)) then
                        jsect = string(1)
                    else
                        jsect = string(x2%p%indx(1)) 
                    endif
                endif

                dim2 = sectors(isect)%ndim
                dim3 = sectors(jsect)%ndim

                if ( .not. x1%p%lsave(isect,this_level,1) ) then
                    call cat_subprod(skiplists,this_level, x1%p, x2%p, csize, string, &
                              index_t_loc, tmp_mat, x1%p%subprod(isect,this_level)%submat(:,:,1))
                    x1%p%lsave(isect, this_level,1) = .true.
                    x1%p%lcopy(isect, this_level) = .true.

                    if (counter == 1) then
                        total_prod = x1%p%subprod(isect,this_level)%submat(:,:,1)
                    else
                        call dgemm( 'N','N', dim3, dim1, dim2, one,                             &
                                    x1%p%subprod(isect,this_level)%submat(:,:,1), max_dim_sect, &
                                    total_prod,                                   max_dim_sect, &
                                    zero, tmp_mat,                                max_dim_sect   )
                        total_prod = tmp_mat
                        num_prod = num_prod + 1
                    endif

                else

                    if (counter == 1) then
                        total_prod = x1%p%subprod(isect,this_level)%submat(:,:,2)
                    else
                        call dgemm( 'N','N', dim3, dim1, dim2, one,                             &
                                    x1%p%subprod(isect,this_level)%submat(:,:,2), max_dim_sect, &
                                    total_prod,                                   max_dim_sect, &
                                    zero, tmp_mat,                                max_dim_sect   )
                        total_prod = tmp_mat
                        num_prod = num_prod + 1
                    endif

                endif

                if (associated(x2%p, tail)) EXIT

                x1%p => x2%p
                x2%p => x2%p%forward(this_level,1)%p

            enddo  
        endif
     
        return
     end subroutine cat_subprod

     subroutine cat_subprod2(skiplists, csize, string, index_t_loc, total_prod)
        type(t_skiplists), pointer :: skiplists
        integer, intent(in)   :: csize
        integer, intent(in)   :: string(csize+1)
        integer, intent(in)   :: index_t_loc(mkink)
        real(dp), intent(inout) :: total_prod(max_dim_sect, max_dim_sect)

        real(dp) :: tmp_mat(max_dim_sect, max_dim_sect)
        type(t_pnode) :: x1, x2, x3, x4
        integer :: this_level
        integer :: counter
        integer :: isect, jsect
        integer :: msect, nsect
        integer :: dim1, dim2, dim3, dim4, dim5, dim6
        integer :: istart
        integer :: vt, vf
        integer :: i,j

        ! the first level
        dim1 = sectors(string(1))%ndim
        x1%p => skiplists%head
        do while(.not. associated(x1%p, skiplists%tail))
            
            if (x1%p%indx(1) == 0) then
                isect = string(1) 
                jsect = string(1) 
            else
                isect = string(x1%p%indx(1))
                jsect = string(x1%p%indx(1)+1)
            endif

            dim2 = sectors(isect)%ndim
            dim3 = sectors(jsect)%ndim

            if (.not. x1%p%lsave(isect,1,1)) then
                x1%p%subprod(isect,1)%submat(:,:,1) = zero
                if (x1%p%indx(1) == 0) then
                    do i=1, dim2
                        x1%p%subprod(isect,1)%submat(i,i,1) = one 
                    enddo
                else
                    istart = sectors(isect)%istart
                    vt = type_v( index_t_loc(x1%p%indx(1)) )
                    vf = flvr_v( index_t_loc(x1%p%indx(1)) ) 
                    do i=1, dim2
                        do j=1, dim3
                            x1%p%subprod(isect,1)%submat(j,i,1) = sectors(isect)%myfmat(vf,vt)%item(j,i) * &
                            expt_v( istart+i-1, index_t_loc(x1%p%indx(1)) )
                        enddo 
                    enddo
                    num_prod = num_prod + 1
                endif
                x1%p%lsave(isect,1,1) = .true.
                x1%p%lcopy(isect,1) = .true.
            endif 
            x1%p => x1%p%forward(1,1)%p
        enddo
 
        do i=2, skiplists%list_level(1)
            x1%p => skiplists%head
            x2%p => skiplists%head%forward(i,1)%p 
            do while(.not. associated(x1%p, skiplists%tail))
                if (x1%p%indx(1) == 0) then
                    isect = string(1) 
                    if (associated(x2%p, skiplists%tail)) then
                        jsect = string(1)
                    else
                        jsect = string(x2%p%indx(1)) 
                    endif
                else
                    isect = string(x1%p%indx(1))
                    if (associated(x2%p, skiplists%tail)) then
                        jsect = string(1)
                    else
                        jsect = string(x2%p%indx(1)) 
                    endif
                endif

                dim2 = sectors(isect)%ndim
                dim3 = sectors(jsect)%ndim

                if ( .not. x1%p%lsave(isect,i,1) ) then
                    counter = 0
                    x3%p => x1%p
                    x4%p => x3%p%forward(i-1,1)%p
                    dim4 = dim2 
                    do while(.not. associated(x3%p, x2%p)) 
                        counter = counter + 1
                        if (x3%p%indx(1) == 0) then
                            msect = string(1) 
                            if (associated(x4%p, skiplists%tail)) then
                                nsect = string(1)
                            else
                                nsect = string(x4%p%indx(1)) 
                            endif
                        else
                            msect = string(x3%p%indx(1))
                            if (associated(x4%p, skiplists%tail)) then
                                nsect = string(1)
                            else
                                nsect = string(x4%p%indx(1)) 
                            endif
                        endif

                        dim5 = sectors(msect)%ndim
                        dim6 = sectors(nsect)%ndim

                        if (counter == 1) then
                            if (x3%p%lcopy(msect, i-1)) then
                                x1%p%subprod(isect,i)%submat(:,:,1) = x3%p%subprod(msect,i-1)%submat(:,:,1)
                            else
                                x1%p%subprod(isect,i)%submat(:,:,1) = x3%p%subprod(msect,i-1)%submat(:,:,2)
                            endif
                        else
                            if (x3%p%lcopy(msect, i-1)) then
                                call dgemm( 'N','N', dim6, dim4, dim5, one,                        &
                                            x3%p%subprod(msect,i-1)%submat(:,:,1),   max_dim_sect, &
                                            x1%p%subprod(isect,i)%submat(:,:,1),     max_dim_sect, &
                                            zero, tmp_mat,                           max_dim_sect   )
                            else
                                call dgemm( 'N','N', dim6, dim4, dim5, one,                        &
                                            x3%p%subprod(msect,i-1)%submat(:,:,2),   max_dim_sect, &
                                            x1%p%subprod(isect,i)%submat(:,:,1),     max_dim_sect, &
                                            zero, tmp_mat,                           max_dim_sect   )

                            endif
                            x1%p%subprod(isect,i)%submat(:,:,1) = tmp_mat
                            num_prod = num_prod + 1
                        endif

                        if (associated(x4%p, x2%p)) EXIT
                        x3%p => x4%p
                        x4%p => x4%p%forward(i-1,1)%p
                    enddo

                    x1%p%lsave(isect,i,1) = .true.
                    x1%p%lcopy(isect,i) = .true.

                endif

                if (associated(x2%p, skiplists%tail)) EXIT
                x1%p => x2%p
                x2%p => x2%p%forward(i,1)%p
            enddo

        enddo

        total_prod = skiplists%head%subprod(string(1), skiplists%list_level(1))%submat(:,:,1)

        return
     end subroutine cat_subprod2

     subroutine random_level(level)
        use spring
        implicit none
 
        integer, intent(out) :: level
        real(dp) :: x

        call random_number(x)
        !x = spring_sfmt_stream()
        x = 2.0_dp + (2**mlevl-2.0_dp) * x 
        
        level = floor(log(x)/log(2.0_dp))
 
        if ( level <=1 ) level = 1
        if ( level > mlevl-1) level = mlevl-1 

        level = mlevl - level 

        return
     end subroutine random_level

  end module m_skiplists
