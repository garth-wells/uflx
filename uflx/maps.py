from abc import ABC, abstractmethod


class AbstractReferenceMap(ABC):
    @abstractmethod
    def push_forward_symbolic(self, function):
        pass

    @abstractmethod
    def pull_back_symbolic(self, function):
        pass


class IdentityReferenceMap(AbstractReferenceMap):
    def push_forward_symbolic(self, function):
        return function

    def pull_back_symbolic(self, function):
        return function
