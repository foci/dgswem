---
trigger: always_on
---

#Fortran Style Suggestions

## Names
- References to variables, subroutines, functions: `snake_case`
- Constants: `CAPS_SNAKE_CASE`
- Types: `PascalCase` 
- Imports use lowercase
- Declaring variables, subroutines, functions, etc.. use an all-caps version of the declaration line, including closing scope.
- Most keywords: `lowercase`

## Comments & documentation
Use Ford/fortls style `!!` under or next to signatures, variables, etc. to document. 

## Imports 
To keep namespaces clear, any imported keyword should be declared during the import statement using `only`

## Alignment 
Indent-to-scope. Break long argument lists with line continuation characters, which are aligned on the right margin. 

##Example
`MODULE EXAMPLE_MODULE
  !! This is an example module 
  use imported_mod, only: var1,           &
                          var2,           &
                          subroutine3
  
  implicit none 
  
  contains

    SUBROUTINE EXAMPLE_SUBROUTINE1(arg)
      !! An example subroutine that calls another
      INTEGER, INTENT(IN) :: ARG 
      INTEGER :: i       

      call example_subroutine2(arg)
      do i = 1, 3
        call subroutine3()
      end do   
    end subroutine EXAMPLE_SUBROUTINE1
    
    SUBROUTINE EXAMPLE_SUBROUTINE2(arg)
      INTEGER, INTENT(IN) :: ARG
    end subroutine EXAMPLE_SUBROUTINE2
end module EXAMPLE_MODULE