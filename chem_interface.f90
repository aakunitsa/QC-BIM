module chem_interface
use qc_system
implicit none

contains

subroutine esp_wrt (im)
    integer, intent(in) :: im
    if (code(1:3) .eq. 'g09') then
        call g09_esp_wrt (im)
    else if (code(1:3) .eq. 'psi') then
        !call psi_esp_wrt (im)
        call qc_nwchem_esp_wrt (im)
    else
        call qc_nwchem_esp_wrt (im)
    end if
end subroutine esp_wrt

subroutine dipole_wrt (im)
    integer, intent(in) :: im
    if (code(1:3) .eq. 'g09') then
        call g09_dipole_wrt (im)
    else if (code(1:3) .eq. 'psi') then
        stop "not implemented"
    else
        call qc_nwchem_dipole_wrt (im)
    end if
end subroutine dipole_wrt

subroutine atomdipole_wrt (im)
    integer, intent(in) :: im
    if (code(1:3) .eq. 'g09') then
        call g09_atomdipole_wrt (im)
    else if (code(1:3) .eq. 'psi') then
        stop "not implemented"
    else
        stop "not yet"
    end if
end subroutine atomdipole_wrt

subroutine esp_read (nsite, chg0)
    integer, intent(in) :: nsite
    real*8, intent(out) :: chg0(nsite)
    if (code(1:3) .eq. 'g09') then
        call g09_esp_read (nsite, chg0)
    else if (code(1:3) .eq. 'psi') then
        !call psi_esp_read (nsite, chg0)
        call qc_nwchem_esp_read (nsite, chg0)
    else
        call qc_nwchem_esp_read (nsite, chg0)
    end if
end subroutine esp_read

subroutine dipole_read (dipole0)
    real*8, intent(out) :: dipole0(3)
    if (code(1:3) .eq. 'g09') then
        call g09_dipole_read (dipole0)
    else if (code(1:3) .eq. 'psi') then
        stop "not implemented"
    else
        call qc_nwchem_dipole_read (dipole0)
    end if
end subroutine dipole_read

subroutine atomdipole_read (nsite, charges, dipoles)
    integer, intent(in) :: nsite
    real*8, intent(out) :: charges(nsite)
    real*8, intent(out) :: dipoles(3, nsite)
    if (code(1:3) .eq. 'g09') then
        call g09_atomdipole_read (nsite, charges, dipoles)
    else if (code(1:3) .eq. 'psi') then
        stop "not implemented"
    else
        stop "not yet"
    end if
end subroutine atomdipole_read

subroutine monomer_wrt (id)
    integer, intent(in) :: id
    if (code(1:3) .eq. 'g09') then
        call g09_monomer_wrt (id)
    else if (code(1:3) .eq. 'psi') then
        call psi_monomer_wrt (id)
    else
        call qc_nwchem_monomer_wrt(id)
    end if
end subroutine monomer_wrt

subroutine monomerbq2_wrt (im, ip, mono, ghost)
    integer, intent(in) :: im, ip, mono
    logical, intent(in) :: ghost
    if (code(1:3) .eq. 'g09') then
        call g09_monomerbq2_wrt (im, ip, mono, ghost)
    else if (code(1:3) .eq. 'psi') then
        call psi_monomerbq2_wrt (im, ip, mono, ghost)
    else
        call qc_nwchem_monomerbq2_wrt (im, ip, mono, ghost)
    end if
end subroutine monomerbq2_wrt

subroutine pot_grad_read (natm, upot, dRR0, dRR1) ! AK
    integer, intent(in)  :: natm
    real*8,  intent(out) :: upot
    real*8,  intent(out) :: dRR0(3,natm)
    real*8, optional, intent(out) :: dRR1(3,natm) ! AK
    if (code(1:3) .eq. 'g09') then
        call g09_pot_grad_read(natm, upot, dRR0)
    else if (code(1:3) .eq. 'psi') then
        if (.not. present(dRR1)) then
            call psi_pot_grad_read(natm, upot, dRR0)
        else
            call psi_pot_grad_read_respa(natm, upot, dRR0, dRR1)
        end if
    else
        call qc_nwchem_pot_grad_read(natm, upot, dRR0)
    end if
end subroutine pot_grad_read

subroutine bq_grad_read (d_bq, im)
    integer, intent(in) :: im
    real*8, intent(out) :: d_bq(3, 2*maxnum_bq*3)
    if (code(1:3) .eq. 'g09') then
        call g09_bq_grad_read(d_bq, im)
    else if (code(1:3) .eq. 'psi') then
        call psi_bq_grad_read(d_bq, im)
    else
        call qc_nwchem_bq_grad_read(d_bq, im)
    end if
end subroutine bq_grad_read

subroutine bq2_grad_read (d_bq)
    real*8, intent(out) :: d_bq(3, 2*maxnum_bq*3)
    if (code(1:3) .eq. 'g09') then
        call g09_bq2_grad_read(d_bq)
    else if (code(1:3) .eq. 'psi') then
        call psi_bq2_grad_read(d_bq)
    else
        call qc_nwchem_bq2_grad_read(d_bq)
    end if
end subroutine bq2_grad_read

subroutine pot_read (upot)
    real*8,  intent(out) :: upot
    if (code(1:3) .eq. 'g09') then
        call g09_pot_read(upot)
    else if (code(1:3) .eq. 'psi') then
        call psi_pot_read(upot)
    else
        call qc_nwchem_pot_read(upot)
    end if
end subroutine pot_read

subroutine dimer_wrt (im, ip)
    integer, intent(in) :: im, ip
    if (code(1:3) .eq. 'g09') then
        call g09_dimer_wrt(im, ip)
    else if (code(1:3) .eq. 'psi') then
        call psi_dimer_wrt(im, ip)
    else
        call qc_nwchem_dimer_wrt(im, ip)
    end if
end subroutine dimer_wrt

subroutine read_hess(nsite, hess0)
    integer, intent(in) :: nsite
    real*8, intent(out) :: hess0(171)
    if (code(1:3) .eq. 'g09') then
        call g09_read_hess (nsite, hess0)
    else if (code(1:3) .eq. 'psi') then
        call psi_read_hess (nsite, hess0)
    else
        call nw_read_hess (nsite, hess0)
    end if
end subroutine read_hess

subroutine read_ddipole(nsite, ddipole0)
    integer, intent(in) :: nsite
    real*8, intent(out) :: ddipole0(54)
    if (code(1:3) .eq. 'g09') then
        call g09_read_ddipole (nsite, ddipole0)
    else if (code(1:3) .eq. 'psi') then
        call psi_read_ddipole (nsite, ddipole0)
    else
        call nw_read_ddipole (nsite, ddipole0)
    end if
end subroutine read_ddipole

end module chem_interface
