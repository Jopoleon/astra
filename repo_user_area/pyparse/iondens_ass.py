class NIAS:

    iondensassign = \
'''! **** Ion density assignment
call markloc("NI assignment")

do J=1, NA1
NI(J) = F1(J) + F2(J) + F3(J) + F4(J) + F5(J) + F6(J) + F7(J) + F8(J) + F9(J)  ! complete AUG
enddo

'''

