import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("wbi", cls=MyCommand)
@click_job_options
@click.option(
    "-v",
    "--version",
    "nbo_version",
    type=click.Choice(["3.1", "7"]),
    default="3.1",
    show_default=True,
    help="NBO version to use for Wiberg bond indices.",
)
@click.pass_context
def wbi(ctx, nbo_version, **kwargs):
    """CLI subcommand for running Gaussian WBI jobs."""

    # get jobrunner for running Gaussian WBI jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    wbi_settings = project_settings.wbi_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    wbi_settings = wbi_settings.merge(job_settings, keywords=keywords)
    if nbo_version == "7" and wbi_settings.route_to_be_written is not None:
        raise click.UsageError(
            "For a custom route, specify pop=nbo7read directly in the route."
        )
    wbi_settings.nbo_version = nbo_version
    check_charge_and_multiplicity(wbi_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(
        f"Weiberg bond index job settings from project: "
        f"{wbi_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.wbi import GaussianWBIJob

    return GaussianWBIJob(
        molecule=molecule,
        settings=wbi_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
