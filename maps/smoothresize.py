import kaa.imlib2 as im
import glob

for pic in glob.glob('/var/www/runat.me/tmp/maps/*/*/*.png'):
    z = im.open(pic)
    if z.size == (512,512):
        b = z.scale((256,256))
        b.save(pic)
