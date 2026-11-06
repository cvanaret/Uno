# Copyright (c) 2021-2026 IRT Saint Exupéry, https://www.irt-saintexupery.com
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""Settings for the Uno funnel SQP algorithm."""

from __future__ import annotations

from unopy.gemseo.settings.base_uno_settings import BaseUnoSettings


class UNO_Funnel_SQP_Settings(BaseUnoSettings):  # noqa: N801
    """Settings for the Uno funnel SQP algorithm."""

    _TARGET_CLASS_NAME = "UNO_Funnel_SQP"
