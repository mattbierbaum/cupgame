import kaa.imlib2 as im
import glob

N = 512

#for pic in glob.glob('/var/www/runat.me/tmp/maps/*/*/*.png'):
for pic in glob.glob('/var/www/runat.me/tmp/maps/energyfar/*.png'):
    z = im.open(pic)
    if z.size == (N,N):
        b = z.scale((N/2,N/2))
        b.save(pic)
