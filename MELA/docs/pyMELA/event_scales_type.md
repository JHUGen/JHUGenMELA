# event_scales_type {#event_scales_type}

This is a very simple class defined in TVar.hh, and has space to record the factorization/renormalization schemes for each event. **You really should have no need to use this independently, it is just here so you can see the output for functions like `Mela::getRenFacScaleMode`**.

## Constructor {#event_scale_constructor}

There is one constructor for this type, and it requires 4 inputs

1. The renormalization scheme (which is a `Mela.EventScaleScheme type` from the @ref scale_enum "EventScaleScheme" enum)
2. The factorization scheme (which is a `Mela.EventScaleScheme type` from the @ref scale_enum "EventScaleScheme" enum)
3. The renormalization scale factor (a double)
4. The factorization scale factor (a double)

## Attributes {#event_scale_attributes}

The type has 4 attributes which directly correspond to the constructor.

### renomalizationScheme {#event_scale_renomalizationScheme}

(No the name is not a typo!) This is a `Mela.EventScaleScheme` type (from the @ref scale_enum "EventScaleScheme" enum) that quantifies the renormalization scheme. Access it via `Mela.event_scales_type.renomalizationScheme`.

### factorizationScheme {#event_scale_factorizationScheme}

This is a `Mela.EventScaleScheme` type (from the @ref scale_enum "EventScaleScheme" enum) that quantifies the the factorization scheme. Access it via `Mela.event_scales_type.factorizationScheme`.

### ren_scale_factor {#event_scale_renscale}

This is a number that quantifies the multiplicative scale factor applied to the renormalization scheme. Access it via `Mela.event_scales_type.ren_scale_factor`.

### fac_scale_factor {#event_scale_facscale}

This is a number that quantifies the multiplicative scale factor applied to the factorization scheme. Access if via `Mela.event_scales_type.fac_scale_factor`.
