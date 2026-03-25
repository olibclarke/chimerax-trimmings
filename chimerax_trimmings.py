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


def _model_id_spec_from_keys(model_ids):
    return "#" + ",".join(".".join(str(part) for part in model_id) for model_id in model_ids)


def _model_spec(models):
    return "#" + ",".join(model.id_string for model in models)


def _atomic_models_by_ids(session, model_ids):
    from chimerax.atomic import AtomicStructure

    wanted = {tuple(model_id) for model_id in model_ids}
    return [
        model
        for model in session.models.list(type=AtomicStructure)
        if _model_id_key(model) in wanted
    ]


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


def _visible_atomic_models(session):
    from chimerax.atomic import AtomicStructure

    return [m for m in session.models.list(type=AtomicStructure) if m.display]


def _next_unused_model_id(session, reserved_ids):
    used_ids = {model.id[0] for model in session.models.list()}
    used_ids.update(reserved_ids)
    model_id = 9000
    while model_id in used_ids:
        model_id += 1
    reserved_ids.add(model_id)
    return model_id


def _close_cycle_model_display_clipboards(session):
    from chimerax.core.commands import run

    state = getattr(session, "_cyclemodeldisplay_state", None)
    if state is None:
        return

    clipboard_ids = state.get("clipboard_ids", ())
    if clipboard_ids:
        run(session, f"close {_model_id_spec_from_keys(clipboard_ids)}")


def _create_cycle_model_display_state(session, models):
    from chimerax.core.commands import run

    reserved_ids = set()
    clipboard_ids = []
    for model in models:
        clipboard_model_id = _next_unused_model_id(session, reserved_ids)
        run(
            session,
            f"combine #{model.id_string} close false modelId {clipboard_model_id} name cyclemodeldisplay_clipboard",
        )
        run(session, f"hide #{clipboard_model_id} models")
        clipboard_ids.append((clipboard_model_id,))

    return {
        "model_ids": tuple(_model_id_key(model) for model in models),
        "clipboard_ids": tuple(clipboard_ids),
        "index": 0,
    }


def _apply_model_display_preset(session, models, preset_index):
    from chimerax.core.commands import run

    spec = _model_spec(models)
    commands = [
        f"graphics silhouettes false; ca_and_sidechains {spec}",
        f"graphics silhouettes false; ca_trace {spec}",
        (
            f"graphics silhouettes false; nucleotides {spec} atoms; ~rib {spec}; ~surf {spec}; "
            f"disp {spec}; ~disp @H*&{spec}; style {spec} stick; style ions sphere; style solvent ball; "
            f"size {spec} ballscale 0.2; size {spec} stickradius 0.07"
        ),
        f"graphics silhouettes false; rib {spec}; ~disp {spec}; ~surf {spec}",
        (
            f"graphics silhouettes false; rib {spec}; ~surf {spec}; ~disp {spec}; disp @CA&protein&{spec}; "
            f"disp @P&nucleic&{spec}; style {spec} stick; disp sidechain&{spec}; "
            f"disp ~backbone&nucleic&{spec}; size {spec} stickradius 0.1; "
            f"size {spec} pseudobondradius 0.1"
        ),
        f"graphics silhouettes false; rib {spec}; ~disp {spec}; ~surf {spec}; rainbow {spec} palette RdYlBu-5",
    ]
    run(session, commands[preset_index])


def cycle_model_display(session):
    from chimerax.core.commands import run

    state = getattr(session, "_cyclemodeldisplay_state", None)
    if state is not None:
        models = _atomic_models_by_ids(session, state["model_ids"])
    else:
        models = _visible_atomic_models(session)
    if not models:
        session.logger.info("No atomic models are currently shown.")
        return

    model_ids = tuple(_model_id_key(model) for model in models)
    if state is None or state["model_ids"] != model_ids:
        _close_cycle_model_display_clipboards(session)
        state = _create_cycle_model_display_state(session, models)

    next_index = (state["index"] + 1) % 7
    if next_index == 0:
        clipboard_spec = _model_id_spec_from_keys(state["clipboard_ids"])
        run(session, f"show {clipboard_spec} models")
        for source_id, target_id in zip(state["clipboard_ids"], state["model_ids"]):
            run(
                session,
                f"mcopy {_model_id_spec_from_keys((source_id,))} to {_model_id_spec_from_keys((target_id,))} settings csv",
            )
        run(session, f"hide {clipboard_spec} models")
    else:
        _apply_model_display_preset(session, models, next_index - 1)

    state["index"] = next_index
    session._cyclemodeldisplay_state = state


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
    register(
        "cyclemodeldisplay",
        CmdDesc(
            synopsis="Cycle shown atomic models through saved, CA-trace, ribbon, and related display presets"
        ),
        cycle_model_display,
        logger=logger,
    )


register_command(session.logger)
