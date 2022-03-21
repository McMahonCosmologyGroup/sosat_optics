import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def beam_convolve(x, y, beam, apert1, apert2,plots):

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

    if plots==1:
        plt.figure(figsize = (14,4))
        plt.subplot(121)
        plt.title("WG Amp")
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.pcolormesh(x,y,10*np.log10(abs(disc)/np.max(abs(disc))),vmin = -25,shading = 'auto')
        plt.colorbar(label = 'dB')
        plt.axis("equal")

        plt.subplot(122)
        plt.title("WG Phase")
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.pcolormesh(x,y,np.arctan2(np.imag(disc),np.real(disc)),shading = 'auto')
        plt.axis("equal")
        plt.xlim(-2,2)
        plt.ylim(-2,2)
        plt.colorbar()
        plt.show()    
        plt.figure(figsize = (14,4))
        plt.subplot(121)
        plt.title("WG FFT Amp")
        plt.pcolormesh(x,y,10*np.log10(abs(disc_fft)/np.max(abs(disc_fft))),vmin = -25,shading = 'auto')
        plt.colorbar(label = 'dB')
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.subplot(122)
        plt.title("WG FFT Phase")
        plt.pcolormesh(x,y,np.arctan2(np.imag(disc_fft),np.real(disc_fft)),shading = 'auto')
        plt.colorbar()
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.show()        
        plt.figure(figsize = (14,4))
        plt.subplot(121)
        plt.title("Input Beam Amp")
        plt.pcolormesh(x,y,10*np.log10(abs(beam)/np.max(abs(beam))),vmin = -25,shading = 'auto')
        plt.colorbar(label = 'dB')
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.subplot(122)
        plt.title("Input Beam Phase")
        plt.pcolormesh(x,y,np.arctan2(np.imag(beam),np.real(beam)),shading = 'auto')
        plt.colorbar()
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.show() 

        plt.figure(figsize = (14,4))
        plt.subplot(121)
        plt.title("Input Beam FFT Amp")
        plt.pcolormesh(x,y,10*np.log10(abs(beam_fft)/np.max(abs(beam_fft))),vmin = -25,shading = 'auto')
        plt.colorbar(label = 'dB')
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.subplot(122)
        plt.title("Input Beam FFT Phase")
        plt.pcolormesh(x,y,np.arctan2(np.imag(beam_fft),np.real(beam_fft)),shading = 'auto')
        plt.colorbar()
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.show()        

        plt.figure(figsize = (14,4))
        plt.subplot(121)
        plt.title("Input Beam FFT x WG FFT Amp")
        plt.pcolormesh(x,y,10*np.log10(abs(beam_conv)/np.max(abs(beam_conv))),vmin = -25,shading = 'auto')
        plt.colorbar(label = 'dB')
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.subplot(122)
        plt.title("Input Beam FFT x WG FFT Phase")
        plt.pcolormesh(x,y,np.arctan2(np.imag(beam_conv),np.real(beam_conv)),shading = 'auto')
        plt.colorbar()
        plt.ylabel("cm")
        plt.xlabel("cm")
        plt.axis("equal")
        plt.show()   

    return x, y, beam_final


def b2a(beam, phase):
    """
    FFT angular space to aperture plane.
    """
    beam_complex = (abs(beam)) * np.exp(phase * np.pi / 180.0 * np.complex(0, 1))
    tmp = np.fft.fftshift(beam_complex)
    tmp = np.fft.ifft2(beam_complex)
    aper_field = np.fft.fftshift(tmp)
    aper_phase = np.arctan2(np.imag(aper_field), np.real(aper_field))
    return aper_field, aper_phase


def a2b(apert, phi):
    """
    FFT aperture plane to angular space.
    """
    apert = (abs(apert)) * np.exp(phi * np.pi / 180.0 * np.complex(0, 1))
    tmp_new = np.fft.fftshift(apert)
    beam_complex = np.fft.fft2(tmp_new)

    beam_complex = np.fft.fftshift(beam_complex)
    phase = np.arctan2(np.imag(beam_complex), np.real(beam_complex))
    return beam_complex, phase


def coords_ang_to_spat(theta_x, theta_y, freq):
    """
    Coordinate transformation from angular to spatial.
    """
    ff_ghz = freq * 1e9
    # Get spatial coordinates
    lam = (3 * 10 ** 8) / ff_ghz
    delta_th = abs(np.max(theta_x) - np.min(theta_x)) / (
        len(theta_x) - 1
    )  # increment in azimuthal angle
    delta_th = abs(np.max(theta_y) - np.min(theta_y)) / (
        len(theta_y) - 1
    )  # increment in azimuthal angle

    x_len = len(theta_x)
    y_len = len(theta_y)

    alpha = lam / delta_th  # increment in x
    beta = lam / delta_th

    delta_x = alpha / x_len  # spatial coordinates conversion
    delta_y = beta / y_len

    x_spat = (
        np.linspace(
            -int((len(theta_x) / 2)), int((len(theta_x) / 2)), int((len(theta_x)))
        )
        * delta_x
        * 1e2
    )
    y_spat = (
        np.linspace(
            -int((len(theta_x) / 2)), int((len(theta_x) / 2)), int((len(theta_x)))
        )
        * delta_y
        * 1e2
    )
    x_spat, y_spat = np.meshgrid(x_spat, y_spat)
    return x_spat, y_spat


def coords_spat_to_ang(x, y, freq):
    """
    Coordinate transformation from spatial to angular.
    """
    ff_ghz = freq * 1e9 # frequency in [Hz]
    # Get spatial coordinates
    lam = (3 * 10 ** 8) / ff_ghz # wavelength in [m]
    
    # Resolution in aperture plane [m]
    delta_x = abs(np.max(x) - np.min(x)) / (len(x) - 1)  # increment in x
    delta_y = abs(np.max(y) - np.min(y)) / (len(y) - 1)  # increment in y

    x_len = len(x)
    y_len = len(y)

    alpha = lam / delta_x
    beta = lam / delta_y

    # Conversion for spatial to angular
    delta_th = alpha / x_len
    delta_ph = beta / y_len

    x_ang = (
        np.linspace(
            -int((len(x) / 2)), int((len(x) / 2)), int((len(x)))
        )
        * delta_th
    )
    y_ang = (
        np.linspace(
            -int((len(x) / 2)), int((len(x) / 2)), int((len(x)))
        )
        * delta_ph
    )
    x_ang, y_ang = np.meshgrid(x_ang, y_ang)
    return x_ang, y_ang
