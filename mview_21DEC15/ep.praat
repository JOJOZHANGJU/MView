# EP  -- compute F0 track and echo to stdout
#	
#	usage:  c:\ProgramFiles\Praat\praatcon.exe ep.praat FNAME.wav
#
# this version packs everything into one line so that matlab DOS response variable will pick it up
# (for some reason first LF terminates response from stdout)

form ExtractPitch
	sentence wav_fname
endform

Read from file... 'wav_fname$'
To Pitch... 0.0 75 600

pitchID = selected("Pitch")
Down to PitchTier
pitchtierID = selected("PitchTier")
num_points = Get number of points

accum$ = ""

for k to num_points
	time = Get time from index... k
	hertz = Get value at index... k
	if hertz = undefined
		new$ = "'time' NaN "
	else
		new$ = "'time' 'hertz' "
	endif
	accum$ = accum$ + new$
endfor

fileappend <stdout> 'accum$'
