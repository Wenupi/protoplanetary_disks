#PLOTS

##################################################

#if(plot_opacity == True):
#    p_plot = np.linspace(2.0,6.0,21)
#    a_plot = np.logspace(np.log10(0.001),np.log10(10.),21)
#
#    EXT_plot = EXT(a_plot,p_plot)
#    ALB_plot = ALB(a_plot,p_plot)
#    
#    plt.close()
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    im = ax.imshow(EXT_plot, cmap='hot', origin='lower', interpolation='gaussian',aspect='auto')
#    cbar = fig.colorbar(im,orientation='vertical')
#    cbar.set_label(r"$\chi_{"+str(wl)+"\mathrm{cm}} [\mathrm{cm}^2/\mathrm{g}_{\mathrm{dust}}]$")
#    plt.xticks( np.linspace(0,len(a_plot)-1,len(a_plot)) , np.round_(np.log10(a_plot),decimals=2))
#    plt.yticks( np.linspace(0,(len(p_plot)-1),len(p_plot)), np.round_(p_plot, decimals=2))
#    plt.xticks(rotation=90)
#    plt.yticks(rotation=0)
#    plt.xlabel('$\log (a_{\mathrm{max}} [\mathrm{cm}])$')
#    plt.ylabel('$p$')
#    plt.savefig('Opacity/EXT_nu'+str(np.round(nu/1.e9,2))+'GHz.pdf', bbox_inches='tight')
#    print ('   - Opacity/EXT_nu'+str(np.round(nu/1.e9,2))+'GHz.pdf saved!')
#
#    plt.close()
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    im = ax.imshow(ALB_plot, cmap='hot', origin='lower', interpolation='gaussian',aspect='auto')
#    cbar = fig.colorbar(im,orientation='vertical')
#    cbar.set_label(r"$\omega_{"+str(wl)+"\mathrm{cm}}$")
#    plt.xticks( np.linspace(0,len(a_plot)-1,len(a_plot)) , np.round_(np.log10(a_plot),decimals=2))
#    plt.yticks( np.linspace(0,(len(p_plot)-1),len(p_plot)), np.round_(p_plot, decimals=2))
#    plt.xticks(rotation=90)
#    plt.yticks(rotation=0)
#    plt.xlabel('$\log (a_{\mathrm{max}} [\mathrm{cm}])$')
#    plt.ylabel('$p$')
#    plt.savefig('Opacity/ALB_nu'+str(np.round(nu/1.e9,2))+'GHz.pdf', bbox_inches='tight')
#    print ('   - Opacity/ALB_nu'+str(np.round(nu/1.e9,2))+'GHz.pdf saved!')

##########################################################    
if(plot_sky == True):
    
    plt.close()
    fig , ax = plt.subplots(nrows=2,ncols=2,figsize=(15,12))
    fig.subplots_adjust(hspace=.15,wspace=.1)
    plt.suptitle('$\lambda = %.2f \ \mathrm{cm}; \ i = %.1f \ \mathrm{deg}$'%(wl,inc*180./np.pi))
    ##########################################################    
    if(intensity_log == True):
        im = ax[0][0].imshow(np.log10(Bright*1.e3), cmap='hot', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[0][0])
        cbar.set_label(r'$\log (I_{\nu})\ \mathrm{[mJy/pix]}$')

        im = ax[0][1].imshow(np.log10(conv_Bright*1.e3), cmap='hot', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[0][1])
        cbar.set_label(r'$\log (I_{\nu})\ \mathrm{[mJy/beam]}$')
        beam_d = patches.Circle((-0.75*length/2.,-0.75*length/2),radius=beam*distance/2.,facecolor='w', edgecolor='k')
        ax[0][1].add_patch(beam_d)
        
    else:
        im = ax[0][0].imshow(Bright*1.e3, cmap='hot', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')        
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[0][0])
        cbar.set_label(r'$I_{\nu} \ \mathrm{[mJy/pix]}$')

        im = ax[0][1].imshow(conv_Bright*1.e3, cmap='hot', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')        
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[0][1])
        cbar.set_label(r'$I_{\nu} \ \mathrm{[mJy/beam]}$')
        beam_d = patches.Circle((-0.75*length/2.,-0.75*length/2),radius=beam*distance/2.,facecolor='w', edgecolor='k')
        ax[0][1].add_patch(beam_d)


    ax[0][0].set_ylabel(r'$y \ \mathrm{[au]}$')
    ax[0][0].set_title(r'$\mathrm{Intensity \ Model}$')
    ax[0][0].set_xticks([])   
    ax[0][1].set_yticks([])
    ax[0][1].set_xticks([])            
    ax[0][1].set_title(r'$\mathrm{Convolved \ intensity}$')

    
    ##########################################################
    if(np.nanmax(op_depth) > 1.):
        Contorno = True
    else:
        Contorno = False
    
    if(opacity_log == True):
        im = ax[1][0].imshow(np.log10(op_depth), cmap='hot_r', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[1][0])
        
        if(Contorno == True):
            con = ax[1][0].contour(np.log10(op_depth), levels=[0.0],colors='k',linestyles='dashed',origin='lower', aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)])
            ax[1][0].plot([],[],'k--',label=r'$\tau_{\chi_{\nu}} = 1$')

    else:
        im = ax[1][0].imshow(op_depth, cmap='hot', origin='lower',aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')
        cbar = fig.colorbar(im,orientation='vertical',ax=ax[1][0])
        if(Contorno == True):
            con = ax[1][0].contour(op_depth, levels=[1.0],colors='k',linestyles='dashed',origin='lower', aspect='auto',
                       extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)])
            
    if (scattering == True):
        cbar.set_label(r'$\log (\tau_{\chi_{\nu}})$')
    else:
        cbar.set_label(r'$\log (\tau_{\kappa_{\nu}})$')


    im = ax[1][1].imshow(conv_TB, cmap='hot', origin='lower',aspect='auto',
                         extent=[np.min(x_array),np.max(x_array),np.min(y_array),np.max(y_array)],interpolation='None')
#    ax[1][1].add_patch(beam_d)
    cbar = fig.colorbar(im,orientation='vertical',ax=ax[1][1])
    cbar.set_label(r'$T_{\mathrm{B}} [\mathrm{K/beam}]$')

        
    ax[1][0].legend(loc=3,framealpha=0.2,fontsize=12)
    ax[1][0].set_ylabel(r'')
    ax[1][1].set_yticks([])   
    ax[1][0].set_title(r'$\mathrm{Optical \ depth \ model}$')

    for j in range(2):                
        ax[1][j].set_xlabel(r'$x \ \mathrm{[au]}$')
    
    plt.savefig('Results/img/Diks_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.pdf', bbox_inches='tight')
    print ('   - Results/img/Diks_wl'+str(np.round(wl,2))+'_inc'+str(np.round(inc*180./np.pi,2))+'deg.pdf saved!')


##########################################################    
if(plot_Temperature == True):
    plt.close()
    fig , ax = plt.subplots(nrows=1,ncols=1,figsize=(14,6))
    im = ax.imshow(np.log10(map_T),cmap='hot', origin='lower',aspect='auto',
                       extent=[x_array[0],x_array[nx-1],z_temp[0],z_temp[len(z_temp)-1]],interpolation='None')
    ax.contour(np.log10(map_T), 4,colors='w',linestyles='dashed',origin='lower', aspect='auto',
                       extent=[x_array[0],x_array[nx-1],z_temp[0],z_temp[len(z_temp)-1]])
    cbar = fig.colorbar(im,orientation='vertical',ax=ax)
    cbar.set_label(r'$\log (T_{\mathrm{d}} [\mathrm{K}])$')
    ax.set_xlim(xmin=0.,xmax=Rout)
    ax.set_xlabel(r'$\varpi \ \mathrm{[au]}$',fontsize=16)
    ax.set_ylabel(r'$z \ \mathrm{[au]}$',fontsize=16)    
    plt.savefig('Results/img/TemperatureStructure_Mdot'+str(Mdot*year/Msun)+'_Tstar'+str(Tstar)+'.pdf', bbox_inches='tight')
    print ('   - Results/img/TemperatureStructure_Mdot'+str(Mdot*year/Msun)+'_Tstar'+str(Tstar)+'.pdf')
