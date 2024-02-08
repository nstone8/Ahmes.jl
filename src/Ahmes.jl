module Ahmes
using Tessen, Unitful

export CompiledGeometry, GWLJob, multijob

maxpoints = 200
#if a number is less than this, we will treat it as zero for the purposes
#of stage movement
zerotol = 1E-12

"""
```
writegwl(io,object[;relorigin])
```
Write a GWL representation of `object` to `io`. if `object <: LocalFrame`
the `relorigin` keyword argument providing the offset from the current
stage position and the origin of the `LocalFrame` coordinate system is required.
"""
function writegwl end

#helper function to do the formatting all in one place
function writegwl(io,p::Vector{<:Number})
    println(io,join(p,"\t"))
end

"""
```julia
writegwl(io,hs,z)
```
`writegwl` needs an additional argument for `HatchedSlice` to define the z coordinate
where the slice should be written.
"""
function writegwl(io,hs::HatchedSlice,z)
    hlines=hs.hatchlines
    rawz=ustrip(u"µm",z)
    curpathlen = 0 #how many points are included on the current call to write
    for path in hlines
        for point in path
            writegwl(io,vcat(point,rawz))
            curpathlen += 1
            if curpathlen >= maxpoints
                #execute the write
                println(io,"write")
                #re-add the last point so we pick up where we left off
                writegwl(io,vcat(point,rawz))
                #reset our counter
                curpathlen = 1
            end
        end
        #execute the write for this path and reset the counter
        println(io,"write")
        curpathlen=0
    end
    #pick up here, hlines should now be organized conveniently
end

function writegwl(io,b::Block{HatchedSlice};
                  relorigin::Vector{<:Unitful.Length})
        #need to move the stage to our origin
    (xmove,ymove,zmove) = ustrip.(u"µm",relorigin)
    (xmove > zerotol) ? println(io, "MoveStageX $xmove") : nothing
    (ymove > zerotol) ? println(io, "MoveStageY $ymove") : nothing
    (zmove > zerotol) ? println(io, "AddZDrivePosition $zmove") : nothing
    #write all the slices in the block
    allz = collect(s[1] for s in slices(b)) |> sort
    for z in allz
        theseslices = filter(slices(b)) do (thisz,thisslice)
            thisz == z
        end
        for s in theseslices
            #do our local rotation
            srot = rotate(s[2],rotation(b))
            writegwl(io,srot,z)
        end
    end
    #return our move
    return relorigin
end

function writegwl(io,sb::SuperBlock{HatchedSlice};
                  relorigin::Vector{<:Unitful.Length})
    curpos = -1*relorigin
    for b in blocks(sb)
        #apply our local rotation
        brot = rotate(b,rotation(sb))
        #origin of brot relative to where we are
        blockorigin = origin(brot) - curpos
        curpos += writegwl(io,brot,relorigin=blockorigin)
    end
    #how far have we traveled
    return curpos + relorigin
end

"""
```julia
CompiledGeometry(filepath,geom...;laserpower,scanspeed)
```
Create a 'pre-compiled' GWL script which will write `geom...` starting at the
current stage position. `laserpower` (should have units of power) and `scanspeed`
(units of velocity) are required to set the writing parameters of the enclosed
geometry. A `Tessen.translate` method is provided if writing at the current stage
position is not desired.
"""
struct CompiledGeometry
    origin::Vector{<:Unitful.Length}
    filepath::String
    displacement::Vector{<:Unitful.Length}
    #enforce that the origin has 3 entries
    function CompiledGeometry(o,f,d)
        @assert length(o) == 3 "origin must consist of 3 entries"
        @assert length(d) == 3 "displacement must consist of 3 entries"
        new(o,f,d)
    end
end

#going to want a Tessen.origin method
Tessen.origin(cg::CompiledGeometry) = cg.origin

function CompiledGeometry(filepath,geom...;
                          laserpower,scanspeed)
    #first write a header to set laser power and scan speed
    rawlaserpower = ustrip(u"mW",laserpower)
    rawscanspeed = ustrip(u"µm/s",scanspeed)
    #variable for tracking stage movements in the local frame
    curpos = [0u"µm",0u"µm",0u"µm"]
    open(filepath,"w") do io
        println(io,"LaserPower $rawlaserpower")
        println(io,"ScanSpeed $rawscanspeed")
        #now write out all the geometry
        for g in geom
            relorigin = origin(g) - curpos
            curpos += writegwl(io,g;relorigin)
        end
    end
    CompiledGeometry([0u"µm",0u"µm",0u"µm"],filepath,curpos)
end

"""
```julia
translate(cg,displacement)
```
`translate` a `CompiledGeometry` object. This has the effect of moving the stage by
`displacement` before running the corresponding .gwl script.
"""
function Tessen.translate(cg::CompiledGeometry,displacement::Vector{<:Unitful.Length})
    #we do this just by changing the origin
    CompiledGeometry(cg.origin + displacement,
                     cg.filepath,
                     cg.displacement)
end

function writegwl(io,cg::CompiledGeometry;
                  relorigin::Vector{<:Unitful.Length})
    #first move to `cg`s origin
    (xmove,ymove,zmove) = ustrip.(u"µm",relorigin)
    (xmove > zerotol) ? println(io, "MoveStageX $xmove") : nothing
    (ymove > zerotol) ? println(io, "MoveStageY $ymove") : nothing
    (zmove > zerotol) ? println(io, "AddZDrivePosition $zmove") : nothing
    #include the compiled .gwl code
    println(io, "include $(cg.filepath)")
    #this will result in a stage movement of relorigin+cg.displacement
    return relorigin+cg.displacement
end

"""
```julia
GwlJob(filepath,cg...;stagespeed,interfacepos)
```
Write a .gwl script to `filepath` consisting of the geometry in
`CompiledGeometry` objects. The `stagespeed` (units of velocity)
and `interfacepos` (units of length) keyword arguments are required.
The `struct` returned by this constructor can be used as an argument
to `multijob`. The script generated by this function will always run at
the current stage position, writing at a specified location can be done
using `multijob`.
"""
struct GWLJob
    filepath::String
    function GWLJob(filepath,cg::CompiledGeometry...;
                    stagespeed,interfacepos)
        open(filepath,"w") do io
            #set the stage velocity
            rawstagespeed = ustrip(u"µm/s",stagespeed)
            println(io,"StageVelocity $rawstagespeed")
            #find interface
            rawip = ustrip(u"µm",interfacepos)
            println(io,"FindInterfaceAt $rawip")
            #write all of cg
            curpos = [0u"µm",0u"µm",0u"µm"]
            for cgi in cg
                relorigin = origin(cgi) - curpos
                curpos += writegwl(io,cgi;relorigin)
            end
        end
        new(filepath)
    end
end

"""
```julia
multijob(filepath,gj...;stagespeed)
```
Create a script to execute multiple `GWLJob` scripts back to back.
`gj...` should be a series of `location => job` pairs specifying where
each job should be written in global coordinates
"""
function multijob(filepath,gj...;stagespeed)
    open(filepath,"w") do io
        #set the stage velocity
        rawstagespeed = ustrip(u"µm/s",stagespeed)
        println(io,"StageVelocity $rawstagespeed")
        for (location,job) in gj
            #move to location in global coordinates
            (lx,ly) = ustrip.(u"µm",location)
            println(io,"GlobalGotoX $lx")
            println(io,"GlobalGotoY $ly")
            println(io,"include $(job.filepath)")
        end
    end
end

end # module Ahmes
