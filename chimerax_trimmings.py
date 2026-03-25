def next_model(session):
    from chimerax.model_panel.tool import model_panel

    mp = model_panel(session, "Model Panel")
    mp._next_model()


def previous_model(session):
    from chimerax.model_panel.tool import model_panel

    mp = model_panel(session, "Model Panel")
    mp._previous_model()


def _model_id_key(model):
    return tuple(model.id)


def _toggle_model_set(session, model_type, state_key, empty_message):
    models = list(session.models.list(type=model_type))
    if not models:
        session.logger.info(empty_message)
        return

    visible = [model for model in models if model.display]
    hidden = [model for model in models if not model.display]

    if visible:
        setattr(session, state_key, [_model_id_key(model) for model in visible])
        for model in visible:
            model.display = False
        return

    saved_ids = getattr(session, state_key, None)
    if saved_ids:
        remaining = {_model_id_key(model): model for model in hidden}
        restored = False
        for model_id in saved_ids:
            model = remaining.get(tuple(model_id))
            if model is not None:
                model.display = True
                restored = True
        if restored:
            return

    for model in hidden:
        model.display = True


def toggle_models(session):
    from chimerax.atomic import AtomicStructure

    _toggle_model_set(
        session,
        AtomicStructure,
        "_togglemodels_hidden_ids",
        "No atomic models are open.",
    )


def toggle_maps(session):
    from chimerax.map import Volume

    _toggle_model_set(
        session,
        Volume,
        "_togglemaps_hidden_ids",
        "No volume maps are open.",
    )


def cycle_lighting(session):
    from chimerax.core.commands import run

    presets = ["simple", "soft", "full", "flat", "gentle"]
    current_index = getattr(session, "_cyclelighting_index", -1)
    next_index = (current_index + 1) % len(presets)
    run(session, f"lighting {presets[next_index]}")
    session._cyclelighting_index = next_index


def register_command(logger):
    from chimerax.core.commands import CmdDesc, register

    register(
        "nextmodel",
        CmdDesc(synopsis="Show the next model in the Model Panel"),
        next_model,
        logger=logger,
    )
    register(
        "prevmodel",
        CmdDesc(synopsis="Show the previous model in the Model Panel"),
        previous_model,
        logger=logger,
    )
    register(
        "togglemodels",
        CmdDesc(synopsis="Hide shown atomic models, then restore the same set on the next use"),
        toggle_models,
        logger=logger,
    )
    register(
        "togglemaps",
        CmdDesc(synopsis="Hide shown volume maps, then restore the same set on the next use"),
        toggle_maps,
        logger=logger,
    )
    register(
        "cyclelighting",
        CmdDesc(synopsis="Cycle through lighting presets"),
        cycle_lighting,
        logger=logger,
    )


register_command(session.logger)
