# Covariant_Derivative

Covariant derivaitve is a pacakge I wrote for general relativty covariant derivation manipuations.

In general relativty, a covariant derivative could take many steps depending on the order of the tensors the derivaitive acts on. In order to make it automatic in computation, we need to set up certain rules.

## Type of Tensors

The tensors types are like

$$G^{ab}{}_{cd}{}^{ef}{}_{hi}$$

where the upper indices and the lower indices need to be handled in different ways. In here we use the labeling like the package xAct, that is, the tensor $G^{ab}{}_{cd}{}^{ef}{}_{hi}$ is denoted as,

```
G[-a,-b,c,d,-e,-f,h,i]
```
The upper indices are labeled with a '-', while the lower indices are labeled with the symbol itself.

