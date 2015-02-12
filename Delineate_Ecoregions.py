#==============================================================================
#
# Name:		Delineate_Ecoregions.py
# Auth:		Alexander Peterson
# Date:		2 Jan. 2015
# Desc:		Python script employing the Arcpy module to extract point
#			coordinates based on ecoregion divisions.
#
#==============================================================================


# Import modules.
import arcpy

# Set workspace.
arcpy.env.workspace = "D:/Dropbox/Workspace/Data/Ecoregions.gdb/"

# Overwrite.
arcpy.env.overwriteOutput = 1

# Set variables.
subregion_lyr = "Subregions"
divisions_field = "DIVISION"

maca_shp = "D:/Dropbox/Workspace/Data/MACA_Coords.shp"

# Make feature layers.
subregion_lyr = arcpy.MakeFeatureLayer_management(subregion_lyr, 'subregion_lyr')
maca_lyr = arcpy.MakeFeatureLayer_management(maca_shp, 'maca_lyr')

# Create list of division names.
subregion_cursor = arcpy.UpdateCursor(subregion_lyr)

# Iterate through rows and create list of names.
# for row in subregion_cursor:
#	print(row.getValue(divisions_field))

# Create a unique list of all field values.
divisions_list = list(set([row[0] for row in arcpy.da.SearchCursor(
	subregion_lyr,(divisions_field))]))

# Iterate over field values, select all corresponding features, and create layers.
prefix = ' "DIVISION" = '
for division in divisions_list:

	# Create expression for selecting features.
	query = prefix + '\'' + division + '\''

	# Select division features.
	arcpy.SelectLayerByAttribute_management(subregion_lyr,"NEW_SELECTION",query)

	# Aggregate subregions by dissolving.
	# arcpy.Dissolve_management(subregion_lyr,new_fc,divisions_field)

	# Select MACA coordinates corresponding to divisions.
	arcpy.SelectLayerByLocation_management(maca_lyr,"WITHIN_A_DISTANCE",
		subregion_lyr,0.042,"NEW_SELECTION")

	# Remove spaces and / characters for name.
	file_name = division.replace(' ','').replace('/','') + '_MACA'

	# Open textfile to store lat/lon values.
	out_file = "D:/Dropbox/Workspace/Data/" + file_name + ".txt"
	print out_file
	maca_file = open(out_file,"w")

	# Build search cursor.
	with arcpy.da.SearchCursor(maca_lyr, ["Field1","Field2"]) as cursor:
		for row in cursor:
			lon = str(row[0])
			lat = str(row[1])
			maca_file.write(lon + "," + lat + "\n")

	# Close file.
	maca_file.close()