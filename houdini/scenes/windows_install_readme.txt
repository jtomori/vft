This is a work in progress guide to get this awesome project working on Windows.

1. Make sure the root folder of the project is named "vft-master" or "vft".

2. Find your "houdini.env" file. An Example location of the file is "C:\Users\Matt\Documents\houdini17.5\houdini.env"
In that file, add the project's houdini path there. Example:
HOUDINI_PATH = $HOUDINI_PATH;C:\Users\Matt\Desktop\vft-master\houdini

Now most examples in the scene vft-master\houdini\scenes\testing_scene_mantra_only.hipnc should work. When you get an error on the HDAs, click "reload" on them.

This guide doesn't cover htoa or OSL yet.
