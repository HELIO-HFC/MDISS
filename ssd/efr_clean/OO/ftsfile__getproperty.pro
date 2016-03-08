;----------------------------------------------------------------------------
; FILE::GETPROPERTY
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;
PRO ftsfile::GetProperty, IMAGE=image, HEADER=header, LOCATION=location

    image = *self.image
    header = *self.header
    location = self.location

END