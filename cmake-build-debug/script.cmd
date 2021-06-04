@echo off

set n=%1

fc stage1.txt test/%n%/stage1.txt
fc stage2.txt test/%n%/stage2.txt
fc stage3.txt test/%n%/stage3.txt
fc z_buffer.txt test/%n%/z_buffer.txt