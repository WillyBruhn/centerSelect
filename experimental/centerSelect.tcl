

set start 483
set end 534
set radius 10
if {$argc > 1} {
    puts "processing $argv argc is $argc"

	set start [lindex $argv 0 ]
	set end [lindex $argv 1 ]
	set radius [lindex $argv 2 ]

} else {
    puts "Too few arguments (found $argc). Need a parameters-file\n"
    
} 

set inputFolder "/home/willy/RedoxChallenges/centerSelect/014"

mol new "$inputFolder/014.pqr" type {pqr} first 0 last -1 step 1 waitfor 1

#mol addfile "$inputFolder/014_pot.dx" type dx first 0 last -1 step 1 waitfor 1 0
#mol new "$inputFolder/014_pot.dx" type dx first 0 last -1 step 1 waitfor 1 0


mol representation Isosurface
mol addrep 1
mol new {/home/willy/RedoxChallenges/centerSelect/014/014_pot.dx} type {dx} first 0 last -1 step 1 waitfor 1 volsets {0 }
animate style Loop



puts "selecting active center from $start to $end with radius $radius"

set activeCenter [atomselect 0 "index $start to $end"]

$activeCenter num

proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z} 
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}


proc points_in_cylinder {} {



}


#atomselect1 num


set geom_c [geom_center $activeCenter]

set second_v {0 0 5}

set surf_c [vecadd $geom_c $second_v]

#puts $surf_c

#puts "geometric center is $geom_c"

## puts "$geom_c[1]"

puts "drawing cylinder from $geom_c to $surf_c"
draw cylinder $geom_c $surf_c radius $radius














