# vmd_checks.tcl -- end-to-end checks of the pdbx plugin inside VMD
#
# Usage (headless):
#   vmd -dispdev text -e tests/vmd_checks.tcl -args /path/to/pdbxplugin.so tests/data [extra.cif ...]
#
# The .gz checks need gzip on the PATH and a plugin built with zlib.
#
# Prints CHECK lines (PASS/FAIL) and exits.

set plugin_so [lindex $argv 0]
set data      [lindex $argv 1]
set extra     [lrange $argv 2 end]
set failures 0

proc check {what cond {detail ""}} {
    global failures
    if {$cond} {
        puts "CHECK PASS $what"
    } else {
        puts "CHECK FAIL $what $detail"
        incr failures
    }
}

plugin dlopen $plugin_so
plugin info "mol file reader" pdbx pinfo
check "VMD uses the new pdbx plugin" [expr {$pinfo(majorversion) >= 2}] \
    "found version $pinfo(majorversion).$pinfo(minorversion) by $pinfo(author)"
plugin info "mol file reader" mmcif minfo
check "the mmcif name is this plugin too" \
    [expr {$minfo(majorversion) >= 2 && $minfo(author) eq $pinfo(author)}] \
    "found version $minfo(majorversion).$minfo(minorversion)"

# 1. unknown element symbols no longer abort the load; type extension picks pdbx
set m [mol new [file join $data nonstandard_elements.cif] waitfor all]
set s [atomselect $m all]
check "nonstandard_elements.cif loads 4 atoms" [expr {[$s num] == 4}]
check "unknown symbols become element X" [expr {[$s get element] eq {X X C X}}] [$s get element]
check "original symbol kept in type" [expr {[$s get type] eq {Cg Cg C Dm}}] [$s get type]
set radii [$s get radius]
set ok 1
foreach r $radii { if {$r <= 0} { set ok 0 } }
check "unknown-element atoms get a nonzero radius" $ok $radii
check "atoms are selectable by the original symbol" [expr {[[atomselect $m "type Cg"] num] == 2}]
$s delete
mol delete $m

# 2. multi-model file -> frames; write all frames and read them back
set m [mol new [file join $data models.cif] type pdbx waitfor all]
check "models.cif: 2 matching models load as 2 frames" [expr {[molinfo $m get numframes] == 2}]
set tmp [file join [pwd] vmd_checks_out.cif]
animate write pdbx $tmp beg 0 end -1 waitfor all $m
set m2 [mol new $tmp type pdbx waitfor all]
check "written file has the same frames" [expr {[molinfo $m2 get numframes] == 2}]
set a [atomselect $m all frame 1]
set b [atomselect $m2 all frame 1]
check "written coordinates match (frame 2)" [expr {[$a get {x y z}] eq [$b get {x y z}]}]
check "written names/resids/chains match" [expr {[$a get {name resname resid chain}] eq [$b get {name resname resid chain}]}]
$a delete; $b delete
mol delete $m; mol delete $m2
file delete $tmp

# 3. explicit bonds with alternate locations
set m [mol new [file join $data altloc_bonds.cif] type pdbx waitfor all]
set sgA [atomselect $m "resid 10 and name SG and altloc A"]
set sg20 [atomselect $m "resid 20 and name SG"]
check "disulfide bonded to altloc A conformer" \
    [expr {[lsearch [lindex [$sgA getbonds] 0] [$sg20 get index]] >= 0}]
set zn [atomselect $m "resname ZN"]
check "metal coordination bond present" [expr {[llength [lindex [$zn getbonds] 0]] >= 1}]
$sgA delete; $sg20 delete; $zn delete
mol delete $m

# 4. IHM spheres-only file -> graphics objects
set m [mol new [file join $data ihm_spheres.cif] type pdbx waitfor all]
check "IHM spheres load as graphics" [expr {[llength [graphics $m list]] == 2}] [graphics $m list]
mol delete $m

# 5. label numbering option
set env(VMDPDBXUSELABELIDS) 1
set m [mol new [file join $data label_vs_auth.cif] type pdbx waitfor all]
check "VMDPDBXUSELABELIDS gives label_seq_id numbering" \
    [expr {[[atomselect $m "resname ALA"] get resid] eq {1 2 3}}]
mol delete $m
unset env(VMDPDBXUSELABELIDS)

# 6. chain selection in the file name, and the mmcif name
set m [mol new "[file join $data chains.cif]?chain=A" type mmcif waitfor all]
set s [atomselect $m all]
check "?chain=A with type mmcif: 5 atoms of chain A" \
    [expr {[$s num] == 5 && [lsort -unique [$s get segname]] eq {A}}] [$s get segname]
check "?chain=A: model 2 is a frame for chain A" [expr {[molinfo $m get numframes] == 2}]
$s delete
mol delete $m
set env(MMCIF_CHAIN) "AA a"
set m [mol new [file join $data chains.cif] waitfor all]
check "MMCIF_CHAIN selects chains AA and a" \
    [expr {[lsort -unique [[atomselect $m all] get segname]] eq {AA a}}]
mol delete $m
unset env(MMCIF_CHAIN)

# 7. .cif.gz: read (VMD guesses the type from the last extension only, so
# "type pdbx" is needed), write, read back
set gz [file join [pwd] vmd_checks_models.cif.gz]
if {[catch {exec gzip -c [file join $data models.cif] > $gz} err]} {
    check "gzip available" 0 $err
} else {
    set m [mol new $gz type pdbx waitfor all]
    check "models.cif.gz loads with type pdbx" \
        [expr {[molinfo $m get numatoms] == 3 && [molinfo $m get numframes] == 2}]
    set out [file join [pwd] vmd_checks_out.cif.gz]
    animate write pdbx $out beg 0 end -1 waitfor all $m
    set fh [open $out rb]; binary scan [read $fh 2] H4 magic; close $fh
    check "animate write pdbx to .cif.gz writes gzip" [expr {$magic eq "1f8b"}] $magic
    set m2 [mol new $out type pdbx waitfor all]
    set a [atomselect $m all frame 1]
    set b [atomselect $m2 all frame 1]
    check "the written .cif.gz reads back" [expr {[molinfo $m2 get numframes] == 2 &&
        [$a get {name x y z}] eq [$b get {name x y z}]}]
    $a delete; $b delete
    mol delete $m; mol delete $m2
    file delete $gz $out
}

# 8. BinaryCIF: picked by the .bcif extension, read like the text file
set bcif [file join [file dirname $data] bcif]
set m [mol new [file join $bcif models.bcif] waitfor all]
set t [mol new [file join $data models.cif] waitfor all]
set a [atomselect $m all frame 1]
set b [atomselect $t all frame 1]
check "models.bcif loads by extension like models.cif" \
    [expr {[molinfo $m get numframes] == 2 &&
           [$a get {name resname resid chain x y z}] eq [$b get {name resname resid chain x y z}]}]
$a delete; $b delete
mol delete $m; mol delete $t
set m [mol new "[file join $bcif chains.bcif]?chain=A" type mmcif waitfor all]
check "BinaryCIF with ?chain=A and type mmcif" [expr {[molinfo $m get numatoms] == 5}]
mol delete $m

# 9. optional extra files: chemical components must come with bond orders
foreach f $extra {
    set m [mol new $f type pdbx waitfor all]
    set s [atomselect $m all]
    set orders [join [$s getbondorders]]
    set hasdouble [expr {[lsearch -real $orders 2.0] >= 0}]
    puts "INFO [file tail $f]: [$s num] atoms, [molinfo $m get numframes] frame(s), double bonds: $hasdouble"
    $s delete
    mol delete $m
}

puts "CHECK SUMMARY failures=$failures"
quit
