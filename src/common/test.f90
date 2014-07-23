
  program test
     !!implicit none

     !!call s_print_exception('sdf','ssss')
     !!call s_print_message('sdf','ssss')
     !!call s_print_error('sdf','ssss')

     !!use linkedlist

     !!implicit none
  
! data is stored in data_t
     !!type :: data_t
     !!    real :: x
     !!end type data_t

! a container for storing data_t pointers
     !!type :: data_ptr
     !!    type(data_t), pointer :: p
     !!end type data_ptr

     !!type(list_t), pointer :: list => null()
     !!type(data_ptr) :: ptr

! allocate a new data element
     !!allocate(ptr%p)
     !!ptr%p%x = 2.7183

! initialize the list with the first data element
     !!call list_init(list, transfer(ptr, list_d))
     !!print *, 'Initializing list with data:', ptr%p, list_count(list)

! allocate a second data element
     !!allocate(ptr%p)
     !!ptr%p%x = 0.5772

! insert the second into the list
     !!call list_insert(list, transfer(ptr, list_d))
     !!print *, 'Inserting node with data:', ptr%p, list_count(list)

! retrieve data from the second node and free memory
     !!ptr = transfer(list_get(list_next(list)), ptr)
     !!print *, 'Second node data:', ptr%p
     !!!deallocate(ptr%p)

! retrieve data from the second node and free memory
     !!ptr = transfer(list_get(list_next(list)), ptr)
     !!print *, 'Second node data:', ptr%p
     !!!deallocate(ptr%p)

! retrieve data from the head node and free memory
     !!ptr = transfer(list_get(list), ptr)
     !!print *, 'Head node data:', ptr%p
     !!!deallocate(ptr%p)

! retrieve data from the head node and free memory
     !!ptr = transfer(list_get(list), ptr)
     !!print *, 'Head node data:', ptr%p
     !!!deallocate(ptr%p)

! allocate a third data element
     !!allocate(ptr%p)
     !!ptr%p%x = 10.00

! put this data to the second node
     !!call list_put(list_next(list), transfer(ptr, list_d))

! retrieve data from the second node and free memory
     !!ptr = transfer(list_get(list_next(list)), ptr)
     !!print *, 'Second node data:', ptr%p
     !!!deallocate(ptr%p)

! retrieve data from the head node and free memory
     !!ptr = transfer(list_get(list), ptr)
     !!print *, 'Head node data:', ptr%p
     !!!deallocate(ptr%p)

! free the list
     !!call list_free(list)

     !!use linkedlist
     !!
     !!implicit none
     !!
     !!type data_t
     !!    character(len=32) :: str_key
     !!    character(len=32) :: str_value
     !!end type data_t
     !!
     !!type(data_t), pointer :: p
     !!type(list_t), pointer :: list => null()
     !!type(list_t), pointer :: curr => null()
     !!integer :: i
     !!
     !!allocate(p)
     !!p%str_key = "key0"
     !!p%str_value = "value0"
     !!call list_init(list, transfer(p, list_d))
     !!print *, list_count(list)
     !!
     !!allocate(p)
     !!p%str_key = "key1"
     !!p%str_value = "value1"
     !!call list_insert(list, transfer(p, list_d))
     !!print *, list_count(list)
     !!
     !!allocate(p)
     !!p%str_key = "key2"
     !!p%str_value = "value2"
     !!call list_insert(list, transfer(p, list_d))
     !!print *, list_count(list)
     !!
     !!allocate(p)
     !!p%str_key = "key3"
     !!p%str_value = "value3"
     !!call list_insert(list, transfer(p, list_d))
     !!print *, list_count(list)
     !!
     !!curr => list
     !!do i=1,list_count(list)-1
     !!    curr => list_next(curr)
     !!    p = transfer( list_get(curr), p )
     !!    print *, p%str_key
     !!enddo
     !!
     !!allocate(p)
     !!p%str_key = "keyX"
     !!p%str_value = "valueX"
     !!call list_put(list_next(list_next(list_next(list))), transfer(p, list_d))
     !!print *, list_count(list)
     !!
     !!curr => list
     !!do i=1,list_count(list)-1
     !!    curr => list_next(curr)
     !!    p = transfer( list_get(curr), p )
     !!    print *, p%str_key
     !!enddo
     !!
     !!call list_free(list)
     !!print *, associated(list%next%next), associated(list%data)
     !!curr => list
     !!do i=1,list_count(list)-1
     !!    curr => list_next(curr)
     !!    !!p = transfer( list_get(curr), p )
     !!    !!print *, p%str_key
     !!enddo
  end program test
