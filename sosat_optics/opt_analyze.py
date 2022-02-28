import numpy as np

def beam_convolve(x, y, beam, apert1, apert2):

    #     Define g
    if apert1 < apert2:
        disc = np.cos(np.pi * y / apert2)
        disc = np.where(abs(y) <= apert2 / 2, disc, 0)
        disc = np.where(abs(x) <= apert1 / 2, disc, 0)

    else:
        disc = np.cos(np.pi * x / apert1)
        disc = np.where(abs(x) <= apert1 / 2, disc, 0)
        disc = np.where(abs(y) <= apert2 / 2, disc, 0)

    tmp = np.fft.fftshift(disc)
    tmp = np.fft.fft2(tmp)
    # Matrix G
    disc_fft = np.fft.fftshift(tmp)

    tmp = np.fft.fftshift(beam)
    tmp = np.fft.fft2(tmp)
    beam_fft = np.fft.fftshift(tmp)

    beam_conv = beam_fft * disc_fft

    tmp = np.fft.ifftshift(beam_conv)
    tmp = np.fft.ifft2(tmp)
    beam_final = np.fft.ifftshift(tmp)

    return x, y, beam_final