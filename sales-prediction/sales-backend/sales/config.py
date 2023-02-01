from pydantic import BaseSettings


class Settings(BaseSettings):
    service_name: str = "sentiment-service"
    root_path: str = ""
    path_prefix: str = ""
    # Explanation configuration
    default_explainer: str = "integrated_gradients"
    default_target: int = 4


settings = Settings()
