#!/bin/bash

bunnykiller=../../build/BunnyKiller

file='cornell.xml'
out='renders/cornell'

$bunnykiller -parse-xml-file $file -film-name $out
