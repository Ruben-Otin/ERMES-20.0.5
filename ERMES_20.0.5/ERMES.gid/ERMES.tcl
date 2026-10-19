########################################################################################
# This file is written in TCL lenguage 
# For more information about TCL look at:  http://www.sunlabs.com/research/tcl/
#
# At least two procs must be in this file:
#
#    InitGIDProject dir - Will be called whenever a project is begun to be used.
#                         where dir is the project's directory
#
#    EndGIDProject      - Will be called whenever a project ends to be used.
# 
########################################################################################

proc ErmesCreateToolbar { {type "DEFAULT INSIDELEFT"} } {
    global MyBitmapsNames MyBitmapsCommands MyBitmapsHelp kratosPriv
    set dir $kratosPriv(dir)

    set MyBitmapsNames(0) "
	        icons/Help.gif \ 
			--- \
			icons/RWTE10.gif \ 
            icons/CoaxTEM.gif \
            icons/RobinCoeff.gif \  
            icons/SourceProp.gif \ 
			--- \
			icons/Materials.gif \
			icons/Plasma.gif \ 
            icons/Dirichlet.gif \
            icons/RobinCond.gif \
            icons/Source.gif \ 
			--- \
			icons/Projection.gif \
			icons/ProblemData.gif \
            icons/Results.gif \ 
			--- \
            icons/Mesh.gif \
            icons/Compute.gif \ 
			--- \
            icons/ERMESNameBar.gif"

    set MyBitmapsCommands(0) [list \
        [list -np- WebPageKratos $kratosPriv(WebH)] \
            [ list "" ] \
	    [list -np- GidOpenMaterials RWPort_TE10] \
		[list -np- GidOpenMaterials CoaxialPort_TEM] \
	    [list -np- GidOpenMaterials Robin_condition_coefficients] \
        [list -np- GidOpenMaterials Current_sources_properties] \
		    [ list "" ] \
        [list -np- GidOpenMaterials IHL_materials] \
        [list -np- GidOpenMaterials Plasma] \
	    [list -np- GidOpenConditions Dirichlet_conditions] \
	    [list -np- GidOpenConditions Robin_conditions] \
	    [list -np- GidOpenConditions Current_sources] \
            [ list "" ] \
        [list -np- GidOpenConditions Field_integrals] \
	    [list -np- GidOpenProblemData Solving_parameters] \
	    [list -np- GidOpenProblemData Results] \
	        [ list "" ] \
	    "Meshing generate" \
	    "Utilities Calculate" \
	        [ list "" ] \
	    [list -np- WebPageKratos $kratosPriv(WebE)]]

    set MyBitmapsHelp(0) { 
	    "ERMES user manual" 
		""\
        "Rectangular waveguide - TE10 properties"\
        "Coaxial waveguide - TEM properties"\
        "Robin condition coefficients"\
        "Current sources properties" 
		""\
	    "IHL materials"\
	    "Plasma"\
		"Dirichlet boundary conditions"\
		"Robin boundary conditions"\
		"Current sources"
		""\
        "Results integrals"\
        "Solving parameters"\
	    "Results" 
		""\
        "Generate mesh"\
        "Calculate" 
		""\
        "ERMES webpage"}
    
    # prefix values:
    #          Pre        Only active in the preprocessor
    #          Post       Only active in the postprocessor
    #          PrePost    Active Always

    set prefix Pre
    set name MyBar
    set kratosPriv(toolbarwin) [CreateOtherBitmaps $name "ERMES toolbar" MyBitmapsNames MyBitmapsCommands MyBitmapsHelp $dir ErmesCreateToolbar $type $prefix]
    AddNewToolbar "ERMES bar" ${prefix}${name}WindowGeom ErmesCreateToolbar
}

proc ErmesDestoyToolbar {} {
    global kratosPriv    
    ReleaseToolbar "ERMES bar"
    rename ErmesCreateToolbar ""    
    catch { destroy $kratosPriv(toolbarwin) }
}

proc ErmesChangeMenus {} {
    global kratosPriv
    GiDMenu::InsertOption "Help" [list [_ "About %s v%s" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)]...] end PREPOST [list Splash $kratosPriv(dir)] "" "" insertafter _
    GiDMenu::InsertOption "Help" [list [_ "%s WebPage" $kratosPriv(ProgramName)]] end PREPOST [list VisitWeb $kratosPriv(WebE)] "" "" insertafter _

    GidChangeDataLabel "Data units" ""
    GidChangeDataLabel "Interval" ""
    GidChangeDataLabel "Conditions" ""
    GidChangeDataLabel "Materials" ""
    GidChangeDataLabel "Problem Data" ""
    GidChangeDataLabel "Local axes" ""
    
    GidAddUserDataOptions "---" 2
    GidAddUserDataOptions "Rectangular waveguide - TE10 properties" "GidOpenMaterials RWPort_TE10" 3
    GidAddUserDataOptions "Coaxial waveguide - TEM properties" "GidOpenMaterials CoaxialPort_TEM" 4
    GidAddUserDataOptions "Generic Robin coefficients" "GidOpenMaterials Robin_condition_coefficients" 5
    GidAddUserDataOptions "Current sources properties" "GidOpenMaterials Current_sources_properties" 6
    GidAddUserDataOptions "---" 7
	GidAddUserDataOptions "IHL Materials" "GidOpenMaterials IHL_materials" 8
	GidAddUserDataOptions "Plasma" "GidOpenMaterials Plasma" 9
    GidAddUserDataOptions "Dirichlet boundary conditions" "GidOpenConditions Dirichlet_conditions" 10
    GidAddUserDataOptions "Robin boundary conditions" "GidOpenConditions Robin_conditions" 11
    GidAddUserDataOptions "Current sources" "GidOpenConditions Current_sources" 12
    GidAddUserDataOptions "---" 13
    GidAddUserDataOptions "Results integrals" "GidOpenConditions Field_integrals" 14
    GidAddUserDataOptions "Solving parameters" "GidOpenProblemData Solving_parameters" 15
    GidAddUserDataOptions "Results" "GidOpenProblemData Results" 16

    after 500 GiDMenu::UpdateMenus
}

proc InitGIDProject { dir } {
    # W $dir
    global kratosPriv
    set kratosPriv(dir) $dir
    set kratosPriv(ProgramName) "ERMES" 
    set kratosPriv(VersionNumber) 20.0.5
    set kratosPriv(Web)  "http://www.gidhome.com"
    set kratosPriv(WebH) "$dir/../Documents/Manual/ERMES_20.0.5_Manual.pdf"
    set kratosPriv(WebE) "https://ruben-otin.blogspot.com/2015/04/ruben-otin-software-ruben-otin-april-19.html"
    if { [GiDVersionCmp 7.0] < 0 } {
        WarnWin [_ "%s v%s is not compatible with GiD version lower than 7.0" $kratosPriv(ProgramName) $kratosPriv(VersionNumber)]
        return
    }
    if { ![GidUtils::IsTkDisabled] } {
        Splash $dir 2501
        ErmesChangeMenus
        ErmesCreateToolbar
        # Register the proc to be automatically called when re-creating all menus (e.g. when doing files new)
        GiD_RegisterPluginAddedMenuProc ErmesChangeMenus
    }
}

proc EndGIDProject {} {
    ErmesDestoyToolbar
    GiD_UnRegisterPluginAddedMenuProc ErmesChangeMenus
}

proc WebPageKratos { dir } {
    global kratosPriv

    VisitWeb $dir
}

proc HelpOnkratos { dir } {
    global kratosPriv

    WarnWin [_ "To obtain help for %s v%s, check the lates news in %s" \
			$kratosPriv(ProgramName) $kratosPriv(VersionNumber) $kratosPriv(WebH)] 
}

#  Program Logo and Version

proc Splash { dir { TimeOut 0 } } {
    global GIDDEFAULT
   
    if { [.gid.central.s disable windows] } { return }

    if { [ winfo exist .splash]} {
	    destroy .splash
	    update
    }

    toplevel .splash
        
    set im [image create photo -file [ file join $dir icons/ERMESLogo.gif]]
    set x [expr [winfo screenwidth .splash]/2-[image width $im]/2]
    set y [expr [winfo screenheight .splash]/2-[image height $im]/2]

    wm geom .splash +$x+$y
    wm transient .splash .gid
    wm overrideredirect .splash 1
    pack [label .splash.l -image $im -relief ridge -bd 2]
    
    bind .splash <1> "destroy .splash"
    bind .splash <KeyPress> "destroy .splash"
    raise .splash .gid
    grab .splash
    focus .splash
    update

    if { $TimeOut > 0 } {
        after $TimeOut "if { [ winfo exist .splash] } { destroy .splash  }"
    }
}