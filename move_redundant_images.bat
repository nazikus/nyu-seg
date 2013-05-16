@ECHO OFF
SET _REDUNDANT=redundant_imgs
SET _LABELS=Label_maps
SET _RANSAC=RANSAC_planes


IF EXIST %_REDUNDANT% GOTO L1
	MD %_REDUNDANT%
	ECHO ^'%_REDUNDANT%^' created.
:L1
IF EXIST %_RANSAC% GOTO L2
	MD %_RANSAC%
	ECHO ^'%_RANSAC%^' created.
:L2
IF EXIST %_LABELS% GOTO L3
	MD %_LABELS%
	ECHO ^'%_LABELS%^' created.
:L3


MOVE *_b_planes_*.png %_REDUNDANT%\
MOVE *_e_support_surface_map_*.png %_REDUNDANT%\
MOVE *_f_plane_map_*.png %_REDUNDANT%\

MOVE *_c_ini_label_map_2*.png %_LABELS%\
MOVE *_c_ini_label_map_3*.png %_LABELS%\
MOVE *_c_ini_label_map_4*.png %_LABELS%\
MOVE *_c_ini_label_map_5*.png %_LABELS%\

MOVE *_a_ransac_*.png %_RANSAC%\

echo Press any key to close this window...
pause > nul