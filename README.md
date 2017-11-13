# Globally Resolved Energy Balance (GREB) Climate Model

GREB is a very simple, globally resolved, and computationally cheap energy balance model, which is capable of simulating the main characteristics of global warming. The model bridges the gap between strongly simplified energy balance models and fully coupled 4-dimensional complex CGCMs. 


## Compile and run

    make greb
    ./greb


## References

Dommenget, D., and J. Floeter 2011: Conceptual Understanding of Climate Change with a Globally Resolved Energy Balance Model. Climate dynamics, 2011, 37, 2143-2165. [link](http://users.monash.edu.au/~dietmard/papers/dommenget.and.floeter.greb.paper.cdym2011.pdf)

[Project website](http://users.monash.edu.au/~dietmard/content/GREB/GREB_model.html)

[Code archive](http://users.monash.edu.au/~dietmard/content/GREB/code_files/greb.web-public.tar.zip) (zip, 65Mb)

The github repo [alex-robinson/greb-ucm](https://github.com/alex-robinson/greb-ucm) contains the model with a few modifications, which have been adopted here. They also have R code for data handling and plotting.


## Notes

### Sea ice

**Dommenget & Floeter (2011), 3.6 Sea ice**: The effect of changes in sea ice cover is only considered in the changes of the effective heat capacity c surf , which changes from the oceans mixed layer values to a 2 m water column over a transition temperature interval, see Fig. 3b. The change in heat capacity goes parallel with the change in the albedo (Fig. 3a). Latent heat releases by freezing and melting are neglected.

In the source code, the relevant lines for sea ice are:

```{fortran}
! Ocean: ice -> albedo/heat capacity linear function of T_surf
  where(z_topo < 0. .and. Tsurf <= To_ice1) a_surf = a_no_ice+da_ice      ! ice
  where(z_topo < 0. .and. Tsurf >= To_ice2) a_surf = a_no_ice             ! no ice
  where(z_topo < 0. .and. Tsurf > To_ice1 .and. Tsurf < To_ice2 ) &
&       a_surf = a_no_ice+da_ice*(1-(Tsurf-To_ice1)/(To_ice2-To_ice1))
```


