# EFM  -- compute formant tracks and echo to stdout (MALE voice)
#	
#	usage:  c:\ProgramFiles\Praat\praatcon.exe ef.praat FNAME.wav
#
# this version packs everything into one line so that matlab DOS response variable will pick it up
# (for some reason first LF terminates response from stdout)

form ExtractFormants
	sentence wav_fname
endform

Read from file... 'wav_fname$'
To Formant (burg)... 0 5 5000 0.025 50
Track... 3 500 1450 2475 3465 4455 1 1 1
Down to FormantTier
Down to TableOfReal... yes no

accum$ = ""

mat$ = selected$("TableOfReal")

nRows = Get number of rows

for r to nRows
	t = TableOfReal_'mat$'[r,1]
	f1 = TableOfReal_'mat$'[r,2]
	f2 = TableOfReal_'mat$'[r,3]
	f3 = TableOfReal_'mat$'[r,4]
	accum$ = accum$ + "'t' 'f1' 'f2' 'f3' "
endfor

fileappend <stdout> 'accum$'
