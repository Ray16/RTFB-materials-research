def execution():
    driver = api.load_module("workspace_standard_polymer_driver", ROOT / "scripts/pipeline/refine_polymer_conformers.py")
    driver.SCRIPT = SCRIPT
    campaign = Path(sys.argv[sys.argv.index("--campaign") + 1]).resolve()
    plan = api.read_plan(campaign)
    stage = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "energy"
    if sys.argv[1] == "launch":
        original_slots = driver.additional_host_slots
        required_gib = plan.get("host_required_bytes", 32 * 2**30) / 2**30 if stage == "energy" else plan.get('native_host_budget_GiB', 256)
        driver.additional_host_slots = lambda available, promises, budget: original_slots(
            available, promises, (required_gib if stage == 'thermal' else max(budget, required_gib)))
        atomic_write_json(campaign / stage / "workspace_host_budget.json", dict(
            required_host_GiB=required_gib, scope="Measured DF resource inventory or declared native-thermal budget"))
    if sys.argv[1] == "worker":
        from redox.workflows.resources import verify_gate_owner, host_capacity
        verify_gate_owner()
        host_capacity(2, minimum_ram_gib=plan.get("host_required_bytes", 32 * 2**30) / 2**30
                      if stage == "energy" else plan.get('native_host_budget_GiB', 256))
        witnesses = controls(plan)
        adapter = api.load_module("workspace_campaign_adapter", ADAPTER)
        old_energy, old_native = driver.calculate_energy, api.calculate_native_G

        def energy(folder, task, pilot=False):
            if pilot:
                raise ValueError("This branch reuses its already checked external energy resource point")
            directory = Path(plan["raw_root"]) / "workspace_execution" / (task["key"] + "_energy")
            with adapter.compact_smd_execution(directory, budget_GiB=plan["execution_workspace_GiB"]):
                result = old_energy(folder, task, pilot=False)
            annotate(folder / "energy/tasks" / task["key"] / "status.json", directory, witnesses, plan)
            return result

        def native(folder, record):
            directory = Path(plan["raw_root"]) / "workspace_execution" / (record["task"]["key"] + "_native")
            with adapter.compact_smd_execution(directory, budget_GiB=plan["execution_workspace_GiB"]):
                result = old_native(folder, record)
            annotate(folder / "thermal/tasks" / record["task"]["key"] / "status.json", directory, witnesses, plan)
            return result

        driver.calculate_energy = energy
        api.calculate_native_G = native
    driver.main()
