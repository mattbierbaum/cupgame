import kaa.imlib2 as im
import glob

#for pic in glob.glob('/var/www/runat.me/tmp/maps/*/*/*.png'):
for pic in glob.glob('/var/www/runat.me/tmp/maps/zoom/*.png'):
    z = im.open(pic)
    if z.size == (1024,1024):
        b = z.scale((512,512))
        b.save(pic)
