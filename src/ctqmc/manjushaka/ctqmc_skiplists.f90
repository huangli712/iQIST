!>>> this module define skip lists 
  module m_skiplists
     use constants
     use control
     use context
     use m_sector

     implicit none

!>>> data structure of skip lists 
     type :: t_skiplists
         integer :: list_level
         integer :: num_node
         type(t_node), pointer :: head
         type(t_node), pointer :: tail
     end type t_skiplists

!>>> data structure of a node
     type :: t_node
         integer :: node_level
         integer :: indx
         type(t_pnode), dimension(:), pointer     :: forward => null()
         logical, dimension(:,:), pointer         :: lsave => null() 
         type(t_subprod), dimension(:,:), pointer :: subprod => null()
     end type t_node

!>>> data structure of a pointer to node
     type :: t_pnode
         type(t_node), pointer :: p
     end type t_pnode

!>>> data structure of a subprod
     type :: t_subprod
         real(dp), dimension(:,:), pointer :: item => null()
     end type t_subprod

!>>> a very large integer
     integer, private, parameter :: max_int = 100000000

!>>> global variables, we use two skip lists
     type(t_skiplists), public, save, pointer :: skiplists_a 
     type(t_skiplists), public, save, pointer :: skiplists_b 

     contains

!>>> create a new skiplists
     subroutine new_skiplists(skiplists)
        implicit none

        type(t_skiplists), pointer :: skiplists
  
        type(t_node), pointer :: head, tail
        integer :: i,j

! allocate memory for skiplists
        allocate(skiplists)

! add a head node with maximum level, and a tail node with minimum level
        call new_node(mlevl, head)
        call new_node(1,     tail)

! set head node
        head%indx = 0
        do i=1, mlevel
            head%forward(i)%p => tail
        enddo
        head%lsave = .false. 

! set tail node
        tail%indx = max_int
        tail%forward(1)%p => null()
 
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
        integer :: i

        if (associated(skiplists)) then
            curr%p => skiplists%head
            do while(associated(curr%p))
                next%p => curr%p%forward(1)%p
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

        allocate(node)

        node%node_level = level
        if (level >= 1) then
            allocate( forward(level) )
            allocate( lsave(nsectors, level) )
            allocate( subprod(nsectors, level) ) 
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
            if (associated(node%subprod)) then
                do i=1, node%node_level
                    do j=1, nsectors
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

        type(t_subprod), pointer :: subprod

        if (.not. associated(subprod%item)) then
            allocate(subprod%item(max_dim_sect, max_dim_sect))
            subprod%item = zero
        endif 

        return
     end subroutine new_subprod_item

!>>> destroy a subprod
     subroutine destroy_subprod_item(subprod)
        implicit none

        type(t_subprod), pointer :: subprod

        if (associated(subprod%item)) then
            deallocate(subprod%item)
        endif

        return
     end subroutine destroy_subprod_item

!>>> insert a node 
     subroutine insert_node(skiplists, level, indx)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: level
        integer, intent(in) :: indx

        type(t_pnode) :: update(mlevl), x
        type(t_node), pointer :: node
        integer :: i

        x%p => skiplists%head
        do i=skiplists%list_level, 1, -1
            do while (x%p%forward(i)%p%indx <= indx)
                x%p => x%p%forward(i)%p
            enddo 
! destroy unused subprod
            do j=1, nsectors
                if (associated(x%p%subprod(j,i)%item) then
                    call destroy_subprod_item(x%p%subprod(j,i)%item)
                endif
            enddo
            x%p%lsave(:,i) = .false.
            update(i)%p => x%p
        enddo 

        if (level > skiplists%list_level - 1) then
            do i=skiplists%list_level,level
                update(i)%p => skiplists%head
                skiplists%head%lsave(:,i) = .false.
            enddo
            skiplists%list_level = level + 1
            skiplists%head%forward(skiplists%list_level)%p => skiplists%tail
            skiplists%head%lsave(:, skiplists%list_level) = .false.
        endif 

        call new_node(level, node)
        node%indx = indx+1
        
        do i=level, 1, -1
            node%lsave(:,i) = .false.
            x = update(i)
            node%forward(i)%p => x%p%forward(i)%p
            x%p%forward(i)%p => node
        enddo

        x%p => node%forward(1)%p
        do while( .not. associated(x%p, skiplists%tail) )
            x%p%indx = x%p%indx + 1
            x%p => x%p%forward(1)%p
        enddo

        skiplists%num_node = skiplists%num_node + 1
  
        return
     end subroutine insert_node

!>>> remove a node
     subroutine remove_node(skiplists, indx)
        implicit none

        type(t_skiplists), pointer :: skiplists
        integer, intent(in) :: indx

        type(t_pnode) :: update(max_level), x
        type(t_node), pointer :: node
        integer :: i,j

        x%p => skiplists%head
        do i=skiplists%list_level, 1, -1
            do while (x%p%forward(i)%p%indx < indx)
                x%p => x%p%forward(i)%p
            enddo 
            do j=1, nsectors
                if (associated(x%p%subprod(j,i)%item) then
                    call destroy_subprod_item(x%p%subprod(j,i)%item)
                endif
            enddo
            x%p%lsave(:,i) = .false.
            update(i)%p => x%p
        enddo 

        node => x%p%forward(1)%p

        x%p => node%forward(1)%p
        do while(.not. associated(x%p, skiplists%tail))
            x%p%indx = x%p%indx - 1
            x%p => x%p%forward(1)%p
        enddo

        do i=1, node%node_level
            update(i)%p%forward(i)%p => node%forward(i)%p 
        enddo

        call destroy_node(node)

        skiplists%num_node = skiplists%num_node - 1

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
        do i=skiplists%list_level, 1, -1
            do while (x%p%forward(i)%p%indx <= indx)
                x%p => x%p%forward(i)%p
            enddo 
            do j=1, nsectors
                if (associated(x%p%subprod(j,i)%item) then
                    call destroy_subprod_item(x%p%subprod(j,i)%item)
                endif
            enddo
            x%p%lsave(:,i) = .false.
        enddo 

        return
     end subroutine change_one_node

!>>> change all node
     subroutine change_all_node(skiplists)
        implicit none

        type(t_skiplists), pointer :: skiplists

        type(t_pnode) :: x
        integer :: i,j
 
        do i=1, skiplists%list_level
            x%p => skiplists%head
            do while(associated(x%p))
                do j=1, nsectors
                    if (associated(x%p%subprod(j,i)%item) then
                        call destroy_subprod_item(x%p%subprod(j,i)%item)
                    endif
                enddo
                x%p%lsave(:,i) = .false.
                x%p => x%p%forward(i)%p
            enddo
        enddo

        return
     end subroutine change_all_node 

!>>> deep copy a skip lists
     subroutine deep_copy_skiplists(this, other)
        implicit none

        type(t_skiplists), pointer :: this
        type(t_skiplists), pointer :: other

        type(t_pnode) :: x1,x2
        integer :: level
        integer :: counter
        integer :: i,j

        ! step 1: build a new skiplists
        call new_skiplists(this)
        
        ! step 2: insert the same node of other to this one 
        x1%p => other%head%forward(1)%p
        counter = 0
        do while(.not. associated(x1%p, other%tail))
            level = x1%p%node_level
            call insert_node(this, level, counter)               
            counter = counter + 1
            x1%p => x1%p%forward(1)%p
        enddo
 
        ! step 3: copy data from other to this one
        do i=1, this%list_level 
            x1%p => this%head
            x2%p => other%head
            do while(.not. associated(x1%p, this%tail))
                do j=1, nsectors
                    if (x2%p%lsave(j,i)) then
                        x1%p%lsave(j,i) = x2%p%lsave(j,i)
                        call new_subprod_item(x1%p%subprod(j,i))                     
                        x1%p%subprod(j,i)%item = x2%p%subprod(j,i)%item
                    endif                 
                enddo
                x1%p => x1%p%forward(i)%p
                x2%p => x2%p%forward(i)%p
            enddo
        enddo

        return
     end subroutine deep_copy_skiplists

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
     
        level = skiplists%list_level
        head => skiplists%head
        tail => skiplists%tail
     
        isect = string(1)
        dim1 = sectors(isect)%ndim

! if we don't save its result, allocate memory and calculate it
        if (.not. head%lsave(isect, level)) then
            call new_subprod_item(head%subprod(isect, level)) 
            call cat_subprod(level, head, tail, csize, string, index_t_loc, head%subprod(isect,level)%item)
            head%lsave(isect, level) = .true.
        endif
        right_mat = head%subprod(isect, level)%item

! special treatment of the last time-evolution operator
        indx = sectors(string(1))%istart
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
     subroutine cat_subprod(level, head, tail, csize, string, index_t_loc, total_prod)
        implicit none
     
        integer, intent(in)   :: level
        type(t_node), pointer :: head
        type(t_node), pointer :: tail
        integer, intent(in)   :: csize
        integer, intent(in)   :: string(csize+1)
        integer, intent(in)   :: index_t_loc(mkink)
        real(dp), intent(in)  :: expt_t_loc(ncfgs)
        real(dp), intent(inout) :: total_prod(max_dim_sect, max_dim_sect)

        type(t_pnode) :: x1, x2
        real(dp) :: tmp_mat(max_dim_sect, max_dim_sect)
        integer :: this_level
        integer :: counter
        integer :: isect, jsect
        integer :: dim1, dim2, dim3
        integer :: istart
        integer :: vt, vf
     
        this_level = level - 1
        total_prod = zero

        if (head%indx == 0) then
            dim1 = sectors(string(1))%ndim
        else
            dim1 = sectors(string(head%indx))%ndim
        endif
 
        counter = 0 
        x1%p => head
        x2%p => head%forward(this_level)%p

        if (this_level == 1) then
            do while(.not. associated(x1%p, tail))
                counter = counter + 1
                if (x1%p%indx == 0) then
                    isect = string(1) 
                    jsect = string(1) 
                else
                    isect = string(x1%p%indx)
                    jsect = string(x1%p%indx+1)
                endif

                dim2 = sectors(isect)%ndim
                dim3 = sectors(jsect)%ndim

                if (.not. x1%p%lsave(isect,1)) then
                    call new_subprod_item(x1%p%subprod(isect,1))
                    x1%p%subprod(isect,1)%item = zero

                    if (x1%p%indx == 0) then
                        do i=1, dim2
                            x1%p%subprod(isect,1)%item(i,i) = one 
                        enddo
                    else
                        istart = sectors(isect)%istart
                        vt = type_v( index_t_loc(x1%p%indx) )
                        vf = flvr_v( index_t_loc(x1%p%indx) ) 
                        do i=1, dim2
                            do j=1, dim3
                                x1%p%subprod(isect,1)%item(j,i) = sectors(i)%myfmat(vf,vt)%item(j,i) &
                                                    * expt_t_loc(istart+i-1, index_t_loc(x1%p%indx))
                            enddo 
                        enddo
                        num_prod = num_prod + 1
                    endif
                    x1%p%lsave(isect,1) = .true.
                endif

                if (counter == 1) then
                    total_prod = x1%p%subprod(isect,1)%item
                else
                    call dgemm( 'N','N', dim3, dim1, dim2, one,           &
                                x1%p%subprod(isect,1)%item, max_dim_sect, &
                                total_prod,                 max_dim_sect, &
                                zero, tmp_mat,              max_dim_sect   )
                    total_prod = tmp_mat
                    num_prod = num_prod + 1
                endif

                x1%p => x1%p%forward(1)%p
            enddo

        else
            do while(.not. associated(x1%p, tail) )
                counter =  counter + 1
                if (x1%p%indx == 0) then
                    isect = string(1) 
                    jsect = string(x2%p%indx) 
                else
                    isect = string(x1%p%indx)
                    jsect = string(x2%p%indx)
                endif

                dim2 = sectors(isect)%ndim
                dim3 = sectors(jsect)%ndim

                if ( .not. x1%p%lsave(isect,this_level) ) then
                    call new_subprod_item( x1%p%subprod(isect, this_level) ) 
                    call cat_subprod( this_level, x1%p, x2%p, csize, string, &
                              index_t_loc, x1%p%subprod(isect,this_level)%item)
                    x1%p%lsave(isect, this_level) = .true.
                endif

                if (counter == 1) then
                    total_prod = x1%p%subprod(isect,this_level)%item
                else
                    call dgemm( 'N','N', dim3, dim1, dim2, one,           &
                                x1%p%subprod(isect,1)%item, max_dim_sect, &
                                total_prod,                 max_dim_sect, &
                                zero, tmp_mat,              max_dim_sect   )
                    total_prod = tmp_mat
                    num_prod = num_prod + 1
                endif

                if (associated(x2%p, tail)) EXIT

                x1%p => x2%p
                x2%p => x2%p%forward(this_level)%p

            enddo  
        endif
     
        return
     end subroutine cat_subprod

  end module m_skiplists
