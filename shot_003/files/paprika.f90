!+----------------------------------------------------------------------+
!| Predetermined or Atomistic PRImary Kinetics Algoorithms              |
!| -                -         ---     -        -                        |
!| CONTENTS: routines for applying a precalculated database of          |
!|          primary damage information to the cluster dynamics          |
!|          equations.  Includes the ability to interpolate by          |
!|          depth or apply the appropriate polynomial coefficients      |
!|          if provided or requested.                                   |
!+----------------------------------------------------------------------+

module PAPRIKA
use DASH
use CLOVES
implicit none

    private
    
    integer,parameter                :: maxSources=9
    character*4,parameter            :: pCode="I1.1"

    integer                          :: iter
    
    type sourceprofile
    
        integer                      :: nBins=1, &
                                        nTerms=2, &
                                        polyOrder=1, &
                                        nDim=1, &
                                        axisOrder(maxAxis)=(/(iter,iter=1,maxAxis)/)
        character*120                :: partitionFile="partition.txt", &
                                        damageFile="damage.txt"
        doubleprecision              :: depth=1.d2, &
                                        flux=1.d0
        logical                      :: defaultPartition = .TRUE., &
                                        defaultDamage = .TRUE., &
                                        conserveMass = .TRUE., &
                                        specified = .FALSE., &
                                        meshForm = .FALSE.

    end type sourceprofile

    type(sourceprofile),save         :: profile(maxSources)

    integer                          :: unitNum = 937, &
                                        parUnit = 938, &
                                        totalTerms = 0, &
                                        cp = 1, &
                                        io
                                        
    doubleprecision                  :: rateScale = 1.d0, &
                                        siphonRate = 0.d0

    character*60                     :: siphonFrom = "", &
                                        siphonTo = ""



    public set_damage, paprika_input
contains
    
    subroutine set_damage()
    integer                          :: i
        io=0
        write(message(1),'(A)') "+--------------------------------------+"
        write(message(2),'(A)') "| Primary Damage Configuration Log     |"
        write(message(3),'(A)') "+--------------------------------------+"
        call write_message(setupFile); profile(1)%specified=.TRUE.
        do i=1,maxSources
            if(profile(i)%specified) then
                call configure_profile(i)
                totalTerms=totalTerms+profile(i)%nTerms
            end if
        end do
        call set_rates(totalTerms)
    end subroutine set_damage

    subroutine configure_profile(i)
    integer,intent(in)               :: i
        if(profile(i)%defaultDamage) then
            write(message(1),'(A)') "WARNING: Using default damage generation rates"
            call write_message(setupFile)
            profile(i)%nBins=1
        else
            call configure_damage(i)
        end if
        if(profile(i)%defaultPartition) then
            write(message(1),'(A)') "WARNING: Using default cluster size distribution"
            call write_message(setupFile)
        else
            call configure_partition(i)
        end if
    end subroutine configure_profile

    subroutine set_rates(nT)
    integer,intent(in)               :: nT
    doubleprecision                  :: damageRate(nT,nMesh)
    integer                          :: idList(nT,nAxis), &
                                        gList(nT),n,nMax,i
        damageRate=0.d0; idList=0; gList=0
        n=1
        do i=1,maxSources
            if(.NOT.profile(i)%specified) cycle
            nMax=n+profile(i)%nTerms-1        
            call apply_profile(i,idList(n:nMax,:),gList(n:nMax),damageRate(n:nMax,:))
            n=nMax+1
        end do
        call set_sources(nT,idList,gList,damageRate)
        call compute_dose_rate(nT,idList,damageRate)
    end subroutine set_rates

    subroutine apply_profile(i,idList,gList,dOut)
    integer,intent(in)               :: i
    doubleprecision,intent(out)      :: dOut(profile(i)%nTerms,nMesh)
    integer,intent(out)              :: idList(profile(i)%nTerms,nAxis), &
                                        gList(profile(i)%nTerms)
    doubleprecision                  :: energyDeposition(profile(i)%nBins,nMesh), &
                                        partition(profile(i)%nTerms,profile(i)%nBins), &
                                        survival(profile(i)%nBins)
        energyDeposition=0.d0; partition=0.d0; dOut=0.d0;
        idList=0; gList=0
        if(profile(i)%defaultDamage) then
            energyDeposition=atomDensity
        else
            call interpolate_damage(i,energyDeposition)
        end if
        if(profile(i)%defaultPartition) then
            idList(1,1)=-1; idList(2,1)=1
            gList(1)=group_ID('void',idList(1,1:nAxis))
            gList(2)=group_ID('int110',idList(2,1:nAxis))
            partition=1.d0; survival=1.d0
        else
            call load_partition(i,idList,gList,survival,partition)
        end if
        dOut = generate_rates(profile(i)%nBins,profile(i)%nTerms,energyDeposition,&
                    survival,partition,idList(:,profile(i)%axisOrder(1)))
        if(profile(i)%conserveMass) call check_consistency(profile(i)%nTerms,idList,dOut)
    end subroutine apply_profile    

    subroutine configure_damage(i)
    integer,intent(in)               :: i
    character*30                     :: itemName, junk
        unitNum=937
        open(unitNum,file=profile(i)%damageFile)
        write(message(1),'(A,'//pCode//',2A)') "Reading profile #",i," spatial distribution from: ",trim(profile(i)%damageFile)
        call write_message(setupFile)
        do while(io.EQ.0)
            read(unitNum,*,IOSTAT=io) itemName
            backspace(unitNum)
            if(io.NE.0) exit
            if((itemName(1:1).EQ."#").OR.(itemName(1:1).EQ."!")) then
                itemName = "#"
            end if
            select case(itemName)
            case("#"); read(unitNum,*,IOSTAT=io) junk
            case("interpolant_range")
                read(unitNum,*,IOSTAT=io) junk, profile(i)%depth
                write(message(1),'(A,'//pCode//',3A,1pE16.8)') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%depth
            case("damage_bins") 
                read(unitNum,*,IOSTAT=io) junk, profile(i)%nBins
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%nBins)//')') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%nBins
                call write_message(setupFile)
            case("produced_species")
                read(unitNum,*,IOSTAT=io) junk, profile(i)%nTerms
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%nTerms)//')') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%nTerms
                call write_message(setupFile)
            case("interpolation_order")
                read(unitNum,*,IOSTAT=io) junk, profile(i)%polyOrder
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%polyOrder)//')') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%polyOrder
                call write_message(setupFile)
            case("mesh_number")
                read(unitNum,*,IOSTAT=io) junk, profile(i)%polyOrder
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%polyOrder)//')') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%polyOrder
                profile(i)%meshForm=.TRUE.
            case("flux")
                read(unitNum,*,IOSTAT=io) junk, profile(i)%flux
                write(message(1),'(A,'//pCode//',3A,1pE16.8)') &
                "  Profile #",i," distribution option ",trim(itemName)," set to ",profile(i)%flux
            case("array")
                close(unitNum)
                return
            case default
                read(parUnit,*,IOSTAT=io) junk
                write(message(1),'(3A,'//pCode//',A)') &
                "CONFIGURATION: Ignoring unrecognized option ",trim(itemName)," in profile #",i," damage file"
                call write_message(warnFile)
            end select
            call write_message(setupFile)
        end do
        close(unitNum)
    end subroutine configure_damage
    
    subroutine interpolate_damage(i,damage)
    integer,intent(in)               :: i
    doubleprecision,intent(out)      :: damage(profile(i)%nBins,nMesh)
    integer                          :: j
        if(profile(i)%meshForm) then
            damage=simple_interpolation(i)
        else
            damage=polynomial_interpolation(i)
        end if
        damage=damage*profile(i)%flux
        write(message(1),'(A)') "Damage rate table in units of defects per volume time"
        write(message(2),'(A16,'//int_to_char(profile(i)%nBins)//'I16)') "x",(/(j,j=1,profile(i)%nBins)/)
        call write_message(setupFile)
        do j=1,nMesh
            write(message(1),'(1pE16.8,'//int_to_char(profile(i)%nBins)//'(1pE16.8))') mesh(j)%location, damage(:,j)
            call write_message(setupFile)
        end do
        close(unitNum)
    end subroutine interpolate_damage 

    function simple_interpolation(i) result(damage)
    integer,intent(in)               :: i
    doubleprecision                  :: damage(profile(i)%nBins,nMesh), &
                                        nodes(1:profile(i)%polyOrder), &
                                        rates(1:profile(i)%polyOrder), x
    integer                          :: bin,o,j,n
    character*30                     :: itemName
        open(unitNum,file=profile(i)%damageFile); io=0  
        do while(io.EQ.0)
            read(unitNum,*,IOSTAT=io) itemName
            if(io.NE.0) then
                write(message(1),'(2A)') "CRITICAL: Could not parse ",trim(profile(i)%damageFile)
                call write_message(errFile); stop
            end if
            select case(itemName)
            case("array")
                exit
            case default
            end select
        end do
        do j=1,profile(i)%nBins
            do n=1,profile(i)%polyOrder
                read(unitNum,*) bin,nodes(n),rates(n)
                if(bin.NE.j) then
                    write(message(1),'(A,'//pcode//')') "DAMAGE: Bin mismatch in damage interpolation for profile #",i
                    write(message(2),'(A,I3.3,A,I3.3,A,1pE16.8)') "  Expecting bin ",j,", found bin ",bin," at depth ",nodes(n)
                end if
            end do
            do n=1,nMesh
                o=1; x=mesh(n)%location
                if((x.LE.0.d0).OR.(x.GE.nodes(profile(i)%polyOrder))) then
                    damage(j,n)=0.d0; cycle
                end if
                do while(x.GT.nodes(o+1))
                    o=o+1
                    if(o.EQ.profile(i)%polyOrder-1) exit
                end do
                damage(j,n)=rates(o)+(rates(o+1)-rates(o))/(nodes(o+1)-nodes(o))*(x-nodes(o))
            end do
        end do
    end function simple_interpolation

    function polynomial_interpolation(i) result(damage)
    integer,intent(in)               :: i
    doubleprecision                  :: damage(profile(i)%nBins,nMesh), &
                                        coeffs(0:profile(i)%polyOrder), &
                                        phi(0:profile(i)%polyOrder,nMesh), x
    integer                          :: bin,o,j
    character*30                     :: word,itemName
        open(unitNum,file=profile(i)%damageFile); io=0  
        do while(io.EQ.0)
            read(unitNum,*,IOSTAT=io) itemName
            if(io.NE.0) then
                write(message(1),'(2A)') "CRITICAL: Could not parse ",trim(profile(i)%damageFile)
                call write_message(errFile); stop
            end if
            select case(itemName)
            case("array")
                exit
            case default
            end select
        end do
        do j=1,nMesh
            x = mesh(j)%location
            phi(0,j) = 1.d0;  if(profile(i)%polyOrder.GT.0) phi(1,j) = x
            do o=2,profile(i)%polyOrder
                phi(o,j) = x*phi(o-1,j)
            end do
        end do
        io=0;damage=0.d0
        do while(io.EQ.0)
            read(unitNum,*,IOSTAT=io) word
            if(word.EQ."end") then
                exit
            else
                backspace(unitNum)
                read(unitNum,*,IOSTAT=io) bin, coeffs
                do j=1,nMesh
                    if(.NOT.is_surface_node(j)) damage(bin,j)=damage(bin,j)+sum(phi(:,j)*coeffs)
                end do
            end if
        end do
    end function polynomial_interpolation

    subroutine configure_partition(i)
    integer,intent(in)               :: i
    character*30                     :: itemName, junk, axisLabel(maxAxis) 
    integer                          :: j,n,axisNum
        open(parUnit,file=profile(i)%partitionFile); io=0
        write(message(1),'(A,'//pCode//',2A)') "Reading profile #",i," clustering partition from: ",trim(profile(i)%partitionFile)
        call write_message(setupFile)
        do while(io.EQ.0)
            read(parUnit,*,IOSTAT=io) itemName
            backspace(parUnit)
            if(io.NE.0) exit
            if((itemName(1:1).EQ."#").OR.(itemName(1:1).EQ."!")) then
                itemName = "#"
            end if
            select case(itemName)
            case("#"); read(parUnit,*,IOSTAT=io) junk
            case("damage_bins") 
                read(parUnit,*,IOSTAT=io) junk, j
                if(j.NE.profile(i)%nBins) then
                    write(message(1),'(A,'//pCode//')') &
                    "CRITICAL: Inconsistent damage parameterization - bin number mismatch for profile ",i
                    call write_message(errFile); stop
                end if
            case("produced_species")
                read(parUnit,*,IOSTAT=io) junk, profile(i)%nTerms
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%nTerms)//')') &
                "  Profile #",i," partition option ",trim(itemName)," set to ",profile(i)%nTerms
                call write_message(setupFile)
            case("axis_order")
                read(parUnit,*,IOSTAT=io) junk, j, axisLabel(1:j)
                do n=1,j
                    axisNum=axis_id(axisLabel(n))
                    if(axisNum.GT.0) then
                        profile(i)%axisOrder(n)=axisNum
                    else
                        write(message(1),'(A,'//pCode//',3A)') &
                        "PRIMARY: profile #",i," axis order specification ",trim(axisLabel(n))," not recognized"
                        call write_message(errFile); stop 
                    end if
                end do
                profile(i)%nDim=j  
            case("clustering_axes")
                read(parUnit,*,IOSTAT=io) junk, profile(i)%nDim
                write(message(1),'(A,'//pCode//',3A,I'//int_length(profile(i)%nDim)//')') &
                "  Profile #",i," partition option ",trim(itemName)," set to ",profile(i)%nDim
                call write_message(setupFile)
            case("conserves_mass")
                read(parUnit,*,IOSTAT=io) junk, profile(i)%conserveMass
                write(message(1),'(A,'//pCode//',3A,L)') &
                "  Profile #",i," partition option ",trim(itemName)," set to ",profile(i)%conserveMass
                call write_message(setupFile)
            case("array")
                close(parUnit)
                return
            case default
                read(parUnit,*,IOSTAT=io) junk
                write(message(1),'(3A,'//pCode//',A)') &
                "CONFIGURATION: Ignoring unrecognized option ",trim(itemName)," in profile #",i," partition file"
                call write_message(warnFile)
            end select
        end do
        close(parUnit)
    end subroutine configure_partition
    
    subroutine load_partition(i,idList,gList,survival,partition)
    integer,intent(in)               :: i
    integer,intent(inout)            :: idList(profile(i)%nTerms,nAxis), gList(profile(i)%nTerms)
    doubleprecision,intent(inout)    :: partition(profile(i)%nTerms,profile(i)%nBins), &
                                        survival(profile(i)%nBins)
    integer                          :: j,n,cMax,sizes(profile(i)%nDim)
    character*30                     :: word,itemName,group(profile(i)%nTerms)
        cMax=min(nAxis,profile(i)%nDim)
        idList=0; gList=0; survival=0.d0; partition=0.d0;
        open(parUnit,file=profile(i)%partitionFile); io=0
        do while(io.EQ.0)
            read(parUnit,*,IOSTAT=io) itemName
            if(io.NE.0) then
                write(message(1),'(3A)') "CONFIGURATION: Could not parse ",trim(profile(i)%partitionFile)
                call write_message(errFile); stop
            end if
            if(itemName.EQ."array") exit
        end do
        read(parUnit,*) word, survival
        do j=1,profile(i)%nTerms
            read(parUnit,*,IOSTAT=io) sizes, group(j), partition(j,:)
            do n=1,profile(i)%nDim
                idList(j,profile(i)%axisOrder(n))=sizes(n)
            end do
            gList(j)=group_id(group(j),idList(j,:))
        end do
        close(parUnit)
    end subroutine load_partition

    function generate_rates(nBins,nTerms,energyDeposition,survival,partition,sizes)
    integer,intent(in)               :: nBins,nTerms,sizes(nTerms)
    doubleprecision,intent(in)       :: energyDeposition(nBins,nMesh), &
                                        survival(nBins), partition(nTerms,nBins)
    integer                          :: i,j
    doubleprecision                  :: generate_rates(nTerms,nMesh)
        do i=1,nTerms; do j=1,nMesh
            generate_rates(i,j)=rateScale*sum(partition(i,:)*survival*energyDeposition(:,j))/abs(sizes(i))
        end do; end do
    end function generate_rates

    subroutine set_sources(nTerms,idList,gList,damageRate)
    integer,intent(in)               :: nTerms,idList(nTerms,nAxis)
    integer,intent(inout)            :: gList(nTerms)
    doubleprecision,intent(in)       :: damageRate(nTerms,nMesh)
    integer                          :: i, j, g0, nTotal
    character*2                      :: q
        nTotal=nTerms
        do i=1,nTerms
            if(.NOT.cluster_exists(gList(i),idList(i,:))) then
                do j=1,group_number()
                    if(cluster_exists(j,idList(i,:))) then
                        write(message(1),'(A,I0.5,4A)') "WARNING: Aliasing a cluster at ",idList(i,1)," from ",&
                                        trim(group_name(gList(i)))," to ",trim(group_name(j))
                        call write_message(setupFile)
                        gList(i)=j
                        exit
                    end if
                end do
            end if
            if(group_name(gList(i)).EQ.siphonFrom) nTotal=nTotal+1
        end do
        call set_source_number(nTotal)
        do i=1,nTerms; if(group_name(gList(i)).EQ.siphonFrom) then
            write(message(1),'(A,I0.5,4A)') "Splitting cluster size ",idList(i,1)," source between ",&
                                        trim(group_name(gList(i)))," and ",trim(siphonTo)
            call write_message(setupFile)
        end if; end do
        write(message(1),'(A)') "Cluster equations with source terms"
        call write_message(setupFile)
        q=int_length(maxval(abs(idList)),1)
        do i=1,nTerms
            if(group_name(gList(i)).EQ.siphonFrom) then
                call add_source_term(idList(i,:),gList(i),damageRate(i,:)*(1.d0-siphonRate))
                write(message(1),'(A,I'//int_length(gList(i))//',A,'//int_to_char(nAxis)//'I'//q//')') &
                    "    group ",gList(i),", composition ",idList(i,:)
                g0=group_id(siphonTo,idList(i,:))
                call add_source_term(idList(i,:),g0,damageRate(i,:)*siphonRate)
                write(message(2),'(A,I'//int_length(g0)//',A,'//int_to_char(nAxis)//'I'//q//')') &
                    "    group ",g0,", composition ",idList(i,:)
            else
                call add_source_term(idList(i,1:nAxis),gList(i),damageRate(i,1:nMesh))
                write(message(1),'(A,I'//int_length(gList(i))//',A,'//int_to_char(nAxis)//'I'//q//')') &
                    "    group ",gList(i),", composition ",idList(i,:)
            end if
            call write_message(setupFile)
        end do
    end subroutine set_sources
    
    subroutine check_consistency(nTerms,idList,rates)
    doubleprecision,intent(inout)    :: rates(nTerms,nMesh)
    integer,intent(in)               :: nTerms,idList(nTerms,nAxis)
    doubleprecision                  :: vacRate(nMesh), intRate(nMesh), delta, delta2
    integer                          :: i,x,mass
        intRate=0.d0; vacRate=0.d0
        do i=1,nTerms
            mass=idList(i,1)
            do x=1,nMesh
                if(mass.LT.0) vacRate(x)=vacRate(x)+rates(i,x)*abs(mass)
                if(mass.GT.0) intRate(x)=intRate(x)+rates(i,x)*mass
            end do
        end do
        delta=(trap_rule(intRate)-trap_rule(vacRate))
        if(trap_rule(vacRate).GT.0) delta=delta/(trap_rule(vacRate))
        if(abs(delta).GT.(1.d-10)) then
            do i=1,nTerms
                if(idList(i,1).GT.0) rates(i,:)=rates(i,:)/(1.d0+delta)
            end do
            delta2=(trap_rule(intRate)-trap_rule(vacRate))/trap_rule(vacRate)
            write(message(1),'(A)') "WARNING: Primary damage rates have been modified to conserve mass"
            write(message(2),'(A,1pE16.8,A,1pE16.8)') "    mass leak has changed from ",delta," to ",delta2
            call write_message(setupFile)
        end if
        do x=1,nMesh
            if((intRate(x)+vacRate(x)).GT.epsilon(vacRate(nMesh))) then
                delta=(intRate(x)-vacRate(x))/max(intRate(x),vacRate(x))
            else
                delta = 0.d0
            end if
            if(delta.GT.(1.d-6)) then
                write(message(1),'(A,1pE16.8,A,I'//int_length(x)//')') "INTERNAL: Mass conservation violated by ",&
                        delta," at node ",x
                call write_message(errFile); stop
            end if
        end do
    end subroutine check_consistency

    subroutine compute_dose_rate(nTerms,idList,rates)
    doubleprecision,intent(inout)    :: rates(nTerms,nMesh)
    integer,intent(in)               :: nTerms,idList(nTerms,nAxis)
    doubleprecision                  :: vacRate(nMesh,nAxis), intRate(nMesh,nAxis), delta
    integer                          :: i,j,x,mass
    character*38                     :: words(nAxis)
        intRate=0.d0; vacRate=0.d0
        do i=1,nTerms; do j=1,nAxis
            mass=idList(i,j)
            do x=1,nMesh
                if(mass.LT.0) vacRate(x,j)=vacRate(x,j)+rates(i,x)*abs(mass)
                if(mass.GT.0) intRate(x,j)=intRate(x,j)+rates(i,x)*mass
            end do
        end do; end do
        delta = (trap_rule(intRate(:,1))-trap_rule(vacRate(:,1)))/max(trap_rule(intRate(:,1)),1.d-15)
        write(message(1),'(A)') "Dose Rate Information"
        write(message(2),'(A,1pE16.8)') "    Average DPA rate ", trap_rule(vacRate(:,1))/atomDensity/foilThickness
        write(message(3),'(A,1pE16.8)') "    Maximum DPA rate ", maxval(vacRate(:,1))/atomDensity
        write(message(4),'(A,1pE16.8)') "    Leaked fraction  ", delta
        call write_message(setupFile)       
        write(message(1),'(A)') "Detailed Source Information"
        write(message(2),'(A19,'//int_to_char(nAxis)//'A38)')  "+-----------------+",&
                            ("-------------------------------------+",i=1,nAxis)
        do i=1,nAxis
            write(words(i),'(A26,A11,A1)') trim(axis_name(i))," ","|"
        end do
        write(message(3),'(A,A18,'//int_to_char(nAxis)//'A38)') "|","depth      |",words
        write(message(4),'(A19,'//int_to_char(nAxis)//'A38)') "+-----------------+",&
                            ("-------------------------------------+",i=1,nAxis)
        call write_message(setupFile)
        do x=1,nMesh
            do i=1,nAxis
                write(words(i),"(A1,1pE16.8,A2,1pE16.8)") "|",vacRate(x,i)/atomDensity," ",intRate(x,i)/atomDensity
            end do            
            write(message(1),'(A,1pE16.8,A,'//int_to_char(nAxis)//'A38,A)') "|",mesh(x)%location," ",words,"|"
            call write_message(setupFile)
        end do
        write(message(1),'(A19,'//int_to_char(nAxis)//'A38)') "+-----------------+",&
                            ("-------------------------------------+",i=1,nAxis)
        call write_message(setupFile)
    end subroutine compute_dose_rate
    
    logical function paprika_input(unitNum,word)
    integer,intent(in)                       :: unitNum
    character*(*),intent(in)                 :: word
    character*30                             :: junk, setName
    logical                                  :: exists
        select case(word)
        !input works like this
        !case("keyword_for_a_variable")
        !   read(unitNum) junk, associatedVariable, optionalAdditionalVairables, ...
        case("dose_rate_scale")
            read(unitNum,*) junk, rateScale
        case("primary_siphon")
            read(unitNum,*) junk, siphonFrom, siphonTo, siphonRate
        case("damage_file")
            read(unitNum,*,iostat=io) junk, setName
            profile(cp)%damageFile=connect_file(setName,exists)
            if(.NOT.exists) then
                write(message(1),'(2A)') "CONFIGURATION: Failed to locate damage file ",trim(setName)
                call write_message(errFile); stop
            end if
            profile(cp)%defaultDamage=.FALSE.
        case("partition_file")
            read(unitNum,*,iostat=io) junk, setName
            profile(cp)%partitionFile=connect_file(setName,exists)
            if(.NOT.exists) then
                write(message(1),'(2A)') "CONFIGURATION: Failed to locate partition file ",trim(setName)
                call write_message(errFile); stop
            end if
            profile(cp)%defaultPartition=.FALSE.
        case("add_damage_profile")
            read(unitNum,*,iostat=io) junk, cp
            profile(cp)%specified=.TRUE.
        case("change_flux")
            read(unitNum,*,iostat=io) junk, cp, profile(cp)%flux
        case default
            paprika_input=.FALSE.
            return
        end select
        paprika_input=.TRUE.
    end function paprika_input

end module PAPRIKA
