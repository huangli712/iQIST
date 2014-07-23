
  program test
     use parser
     use linkedlist
 
     !!call p_create("solver.ctqmc.in")
     !!call p_parse()

     type(T_data) :: data
     type(T_node), pointer :: list

     data%str_key = "key_0"
     data%str_value = "value_0"
     call list_create(list, data)
     call list_navigator(list)

     data%str_key = "key_1"
     data%str_value = "value_1"
     call list_insert_head(list, data)
     call list_navigator(list)
     
     data%is_valid = .TRUE.
     data%str_key = "key_2"
     data%str_value = "value_2"
     call list_insert_head(list, data)
     call list_navigator(list)

     !call list_delete_head(list)
     !call list_navigator(list)
     !call list_delete_head(list)
     !call list_navigator(list)
     !call list_delete_head(list)
     !call list_navigator(list)

     !!print *, 'count', list_count(list)
     call list_destroy(list)
     !!print *, associated(list), associated(list%next)
     !!call list_display(1, list%next)
     !!print *, list_count(list)
  end program test
