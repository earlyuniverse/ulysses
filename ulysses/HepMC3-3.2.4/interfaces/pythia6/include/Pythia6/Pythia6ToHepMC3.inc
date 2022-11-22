CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Copyright (C) 2014-2020 The HepMC collaboration
C...These variables will be used for conversion. Block size is different
C...in Pythia6 and in HepMC3, so the most simple portable way is to have
C... a second block of same size as in HepMC3 and  copy the content of
C...block directly.
      INTEGER NMXHEPL
      PARAMETER (NMXHEPL=10000)
      COMMON /HEPEVTL/  NEVHEPL,NHEPL,ISTHEPL(NMXHEPL),IDHEPL(NMXHEPL),
     &           JMOHEPL(2,NMXHEPL),JDAHEPL(2,NMXHEPL),PHEPL(5,NMXHEPL),
     &                 VHEPL(4,NMXHEPL)
      INTEGER          NEVHEPL,NHEPL,ISTHEPL,IDHEPL,JMOHEPL,JDAHEPL
      DOUBLE PRECISION PHEPL,VHEPL


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Glue functions
      INTEGER HepMC3_new_writer,HepMC3_delete_writer
      INTEGER HepMC3_set_hepevt_address
      INTEGER HepMC3_convert_event,HepMC3_write_event,HepMC3_clear_event
      INTEGER HepMC3_set_attribute_int,HepMC3_set_attribute_double
      INTEGER HepMC3_set_pdf_info,HepMC3_set_cross_section
      INTEGER HepMC3_new_weight
      INTEGER HepMC3_set_weight_by_index
      INTEGER HepMC3_set_weight_by_name
      EXTERNAL HepMC3_new_writer,HepMC3_delete_writer
      EXTERNAL HepMC3_set_hepevt_address
      EXTERNAL HepMC3_convert_event,HepMC3_write_event
      EXTERNAL HepMC3_clear_event
      EXTERNAL HepMC3_set_attribute_int,HepMC3_set_attribute_double
      EXTERNAL HepMC3_set_pdf_info,HepMC3_set_cross_section
      EXTERNAL HepMC3_new_weight
      EXTERNAL HepMC3_set_weight_by_index
      EXTERNAL HepMC3_set_weight_by_name
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
