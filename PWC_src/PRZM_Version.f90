
Module PRZM_VERSION
    IMPLICIT NONE

CONTAINS

    SUBROUTINE przm_id

    USE  constants_and_variables, ONLY: Version_Number
    IMPLICIT NONE
    WRITE(*,'(A19,1X,F7.3)')' PRZM5-VVWM Version:', Version_Number
    WRITE(*, *)
    WRITE(*, *) '    For technical support contact:'
    WRITE(*, *) '           Dirk F. Young '
    WRITE(*, *) '      Office of Pesticide Programs '
    WRITE(*, *) 'United States Environmental Protection Agency'
    WRITE(*, *) '      Washington, DC 20460-0001'
    WRITE(*, *) '      E-mail: young.dirk@epa.gov'
    WRITE(*, *)
    END SUBROUTINE przm_id

END MODULE PRZM_VERSION
