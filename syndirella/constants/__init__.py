#!/usr/bin/env python3
"""
Constants module containing enumerators and constants used throughout syndirella.
"""

from enum import Enum, auto


class DatabaseSearchTool(Enum):
    """Enumerator for database search tools."""
    MANIFOLD = "manifold"
    ARTHOR = "arthor"
    # HIPPO = "hippo"  # HIPPO dependency removed
    
    @classmethod
    def from_string(cls, value: str) -> 'DatabaseSearchTool':
        """Create DatabaseSearchTool from string value."""
        try:
            return cls(value.lower())
        except ValueError:
            raise ValueError(f"Invalid database search tool: {value}. Valid options: {[tool.value for tool in cls]}")
    
    def __str__(self) -> str:
        return self.value


class RetrosynthesisTool(Enum):
    """Enumerator for retrosynthesis tools."""
    MANIFOLD = "manifold"
    AIZYNTHFINDER = "aizynthfinder"
    
    @classmethod
    def from_string(cls, value: str) -> 'RetrosynthesisTool':
        """Create RetrosynthesisTool from string value."""
        try:
            return cls(value.lower())
        except ValueError:
            raise ValueError(f"Invalid retrosynthesis tool: {value}. Valid options: {[tool.value for tool in cls]}")
    
    def __str__(self) -> str:
        return self.value


# Default values
DEFAULT_DATABASE_SEARCH_TOOL = DatabaseSearchTool.ARTHOR
DEFAULT_RETROSYNTHESIS_TOOL = RetrosynthesisTool.AIZYNTHFINDER